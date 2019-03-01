#include <fstream>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdarg.h>

#include <cstdio>
#include <cstdlib>

//#include <regex.h>
#include <fnmatch.h>
#include <tk.h>

#include "andrea.h"
#include "particles_to_slices.h"
#include "giovanni.h"
#include "beam.h"
#include "beamline.h"
#include "girder.h"

#include "placet.h"
#include "placet_tk.h"

#include "parallel_tracking.hh"
#include "file_stream.hh"

#include "kalman.hh"
#include "latex_picture.h"

#include "defaults.h"
#include "config.h"

#include "dl_element.h"
#include "crab.h"
#include "link.h"
#include "rmatrix.h"
#include "manual.h"
#include "aml.h"

#include "placeti3.h"

#include "quadbpm.h"
#include "multipole.h"
#include "sbend.h"

#include "placet_interface.hh"

#define EMITT_MINIMIZATION(axis)					\
  double emitt_ ## axis ## _disp_free(const BEAM * beam )		\
  {									\
    int i, j, nb, ns, k;						\
    double sum11, sum12, sum22;                                         \
      								        \
    nb = beam->bunches;						        \
    ns = beam->slices_per_bunch * beam->macroparticles;		        \
    double axis ## m=0.0;						\
    double axis ## pm=0.0;						\
    double axis ## e=0.0;                                               \
    double axis ## pe=0.0;                                              \
    double sume = 0.0;                                                  \
    double ee = 0.0;                                                    \
    double wgtsum = 0.0;						\
      							  	        \
    if (beam->particle_beam) {				    	        \
      ns = beam->slices_per_bunch;					\
      k = 0;								\
      for (i = 0; i < ns * (nb - 1); i++) {				\
	for (j = 0; j < beam->particle_number[i]; ++j) {		\
	  double wgt = fabs(beam->particle[k].wgt);		        \
          double energy = fabs(beam->particle[k].energy);               \
	  wgtsum += wgt;						\
	  axis ## m += wgt * beam->particle[k].axis;			\
	  axis ## pm += wgt * beam->particle[k].axis ## p;		\
          axis ## e += wgt * (beam->particle[k].axis * energy);         \
          axis ## pe += wgt * (beam->particle[k].axis ## p * energy);   \
          ee += wgt * energy * energy;                                  \
	  sume += wgt * energy;		                                \
	  k++;							        \
	}								\
      }								        \
      for (; i < ns * nb; i++) {					\
	for (j = 0; j < beam->particle_number[i]; ++j) {		\
	  double wgt = fabs(beam->particle[k].wgt) * beam->last_wgt;    \
          double energy = fabs(beam->particle[k].energy);               \
	  wgtsum += wgt;						\
	  axis ## m += wgt * beam->particle[k].axis;			\
	  axis ## pm += wgt * beam->particle[k].axis ## p;		\
          axis ## e += wgt * (beam->particle[k].axis * energy);         \
          axis ## pe += wgt * (beam->particle[k].axis ## p * energy);   \
          ee += wgt * energy * energy;                                  \
	  sume += wgt * energy;		                                \
	  k++;							        \
	}								\
      }								        \
      sume /= wgtsum;							\
      axis ## m /= wgtsum;						\
      axis ## pm /= wgtsum;						\
      axis ## e /= wgtsum;                                              \
      axis ## pe /= wgtsum;                                             \
      ee /= wgtsum;                                                     \
      double a = (axis ## e-axis ## m*sume)/(ee-sume*sume);             \
      double b = (axis ## pe-axis ## pm*sume)/(ee-sume*sume);           \
      sum11 = 0.0;							\
      sum12 = 0.0;							\
      sum22 = 0.0;							\
      k = 0;								\
      for (i = 0; i < ns * (nb - 1); i++) {				\
	for (j = 0; j < beam->particle_number[i]; ++j) {		\
	  double wgt = fabs(beam->particle[k].wgt);		        \
          double energy = fabs(beam->particle[k].energy);               \
	  double axis = (beam->particle[k].axis - axis ## m) - a * (energy - sume); \
	  double axis ## p = (beam->particle[k].axis ## p - axis ## pm) - b * (energy - sume); \
	  sum11 += axis * axis * wgt;					\
	  sum12 += axis * axis ## p * wgt;				\
	  sum22 += axis ## p * axis ## p * wgt;			        \
	  k++;							        \
	}								\
      }								        \
      for (; i < ns * nb; i++) {					\
	for (j = 0; j < beam->particle_number[i]; ++j) {		\
	  double wgt = fabs(beam->particle[k].wgt) * beam->last_wgt;    \
          double energy = fabs(beam->particle[k].energy);               \
	  double axis = (beam->particle[k].axis - axis ## m) - a * (energy - sume); \
	  double axis ## p = (beam->particle[k].axis ## p - axis ## pm) - b * (energy - sume); \
	  sum11 += axis * axis * wgt;					\
	  sum12 += axis * axis ## p * wgt;				\
	  sum22 += axis ## p * axis ## p * wgt;			        \
	  k++;							        \
	}								\
      }								        \
    } else {								\
      for (i = 0; i < ns * (nb - 1); i++) {				\
	double wgt = fabs(beam->particle[i].wgt);			\
        double energy = fabs(beam->particle[i].energy);                 \
	wgtsum += wgt;						        \
	axis ## m += wgt * beam->particle[i].axis;			\
	axis ## pm += wgt * beam->particle[i].axis ## p;		\
        axis ## e += wgt * (beam->particle[i].axis * energy);           \
        axis ## pe += wgt * (beam->particle[i].axis ## p * energy);     \
        ee += wgt * energy * energy;                                    \
	sume += wgt * energy;			                        \
      }								        \
      for (; i < ns * nb; i++) {					\
	double wgt = fabs(beam->particle[i].wgt * beam->last_wgt);	\
        double energy = fabs(beam->particle[i].energy);                 \
	wgtsum += wgt;						        \
	axis ## m += wgt * beam->particle[i].axis;			\
	axis ## pm += wgt * beam->particle[i].axis ## p;		\
        axis ## e += wgt * (beam->particle[i].axis * energy);           \
        axis ## pe += wgt * (beam->particle[i].axis ## p * energy);     \
        ee += wgt * energy * energy;                                    \
	sume += wgt * energy;			                        \
      }								        \
      ee /= wgtsum;                                                     \
      sume /= wgtsum;							\
      axis ## m /= wgtsum;						\
      axis ## pm /= wgtsum;						\
      axis ## e /= wgtsum;                                              \
      axis ## pe /= wgtsum;                                             \
      double a = (axis ## e-axis ## m*sume)/(ee-sume*sume);             \
      double b = (axis ## pe-axis ## pm*sume)/(ee-sume*sume);           \
      sum11 = 0.0;							\
      sum12 = 0.0;							\
      sum22 = 0.0;							\
      if (placet_switch.first_order) {				        \
	for (i = 0; i < ns * (nb - 1); i++) {				\
	  double wgt = fabs(beam->particle[i].wgt);			\
          double energy = fabs(beam->particle[i].energy);               \
	  double axis = (beam->particle[i].axis - axis ## m) - a * (energy - sume); \
	  double axis ## p = (beam->particle[i].axis ## p - axis ## pm) - b * (energy - sume); \
	  sum11 += axis * axis * wgt;					\
	  sum12 += axis * axis ## p * wgt;				\
	  sum22 += axis ## p * axis ## p * wgt;			        \
	}								\
	for (; i < ns * nb; i++) {					\
	  double wgt = fabs(beam->particle[i].wgt) * beam->last_wgt;	\
          double energy = fabs(beam->particle[i].energy);               \
	  double axis = (beam->particle[i].axis - axis ## m) - a * (energy - sume); \
	  double axis ## p = (beam->particle[i].axis ## p - axis ## pm) - b * (energy - sume); \
	  sum11 += axis * axis * wgt;					\
	  sum12 += axis * axis ## p * wgt;				\
	  sum22 += axis ## p * axis ## p * wgt;			        \
	}								\
      } else {							        \
	for (i = 0; i < ns * (nb - 1); i++) {				\
	  double wgt = fabs(beam->particle[i].wgt);			\
          double energy = fabs(beam->particle[i].energy);               \
	  double axis = (beam->particle[i].axis - axis ## m) - a * (energy - sume); \
	  double axis ## p = (beam->particle[i].axis ## p - axis ## pm) - b * (energy - sume); \
	  sum11 += (axis * axis + beam->sigma_ ## axis ## axis[i].r11) * wgt; \
	  sum12 += (axis * axis ## p + beam->sigma_ ## axis ## axis[i].r12) * wgt; \
	  sum22 += (axis ## p * axis ## p + beam->sigma_ ## axis ## axis[i].r22) * wgt;	\
	}								\
	for (; i < ns * nb; i++) {					\
	  double wgt = fabs(beam->particle[i].wgt) * beam->last_wgt;	\
          double energy = fabs(beam->particle[i].energy);               \
	  double axis = (beam->particle[i].axis - axis ## m) - a * (energy - sume); \
	  double axis ## p = (beam->particle[i].axis ## p - axis ## pm) - b * (energy - sume); \
	  sum11 += (axis * axis + beam->sigma_ ## axis ## axis[i].r11) * wgt; \
	  sum12 += (axis * axis ## p + beam->sigma_ ## axis ## axis[i].r12) * wgt; \
	  sum22 += (axis ## p * axis ## p + beam->sigma_ ## axis ## axis[i].r22) * wgt;	\
	}                                                               \
      }								        \
    }									\
    sum11 /= wgtsum;							\
    sum12 /= wgtsum;							\
    sum22 /= wgtsum;							\
    return sqrt((sum11 * sum22 - sum12 * sum12)) * sume / EMASS * 1e-12 / EMITT_UNIT; \
  }

namespace {
  
  int tk_DispersionFreeEmittance(ClientData /*clientdata*/,Tcl_Interp *interp,int argc,char *argv[])
  {
    int on=-1, off=-1;
    Tk_ArgvInfo table[]={
      {"-on", TK_ARGV_CONSTANT, (char *) 1, (char *) &on, "turn on disp. free emittance calculation"},
      {"-off", TK_ARGV_CONSTANT, (char *) 1, (char *) &off, "turn off disp. free emittance calculation"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,TK_ARGV_NO_LEFTOVERS))!=TCL_OK){
      return error;
    }
    if (on==-1 && off==-1) {
      Tcl_SetResult(interp,"DispersionFreeEmittance: you should select either '-on' or '-off'", TCL_VOLATILE);
      return TCL_ERROR;
    }
    if (on==1 && off==1) {
      Tcl_SetResult(interp,"DispersionFreeEmittance: 'on' and 'off' should not be present at the same time.", TCL_VOLATILE);
      return TCL_ERROR;
    }
    ::Andrea::UseDispersionFreeEmittance = (off == -1);
    return TCL_OK;
  }

  double get_beam_energy(BEAM *beam )
  {
    double wgtsum=0.0,esum=0.0;
    for (int i=0;i < beam->slices; ++i) {
      const PARTICLE &particle = beam->particle[i];
      double wgt = fabs(particle.wgt);
      wgtsum += wgt;
      esum += wgt * fabs(particle.energy);
    }
    return esum / wgtsum;
  }
  
  int tk_SetReferenceEnergy(ClientData /*clientdata*/,Tcl_Interp *interp,int argc,char *argv[])
  {
    if (argc==2) {
      Andrea::ReferenceEnergy = atof(argv[1]);
      return TCL_OK;
    } else {
      beamline_survey_hook_interp = interp;
      char *beamline_name = NULL;
      char *beam_name = NULL;
      char *survey_name = "None";
      int start = 0;
      int end = -1;
      int iter = 3;
      int error;
      Tk_ArgvInfo table[] = {
        {(char*)"-beamline",TK_ARGV_STRING,(char*)NULL,(char*)&beamline_name,
         (char*)"Name of the beamline to be saved"},
        {(char*)"-beam",TK_ARGV_STRING,(char*)NULL,(char*)&beam_name,
         (char*)"Name of the beam to use for the calculation"},
        {(char*)"-survey",TK_ARGV_STRING,(char*)NULL,(char*)&survey_name,
         (char*)"Type of prealignment survey to be used, defaults to None"},
	{(char*)"-start",TK_ARGV_INT,(char*)NULL,(char*)&start,
	 (char*)"First element to be considered"},
	{(char*)"-end",TK_ARGV_INT,(char*)NULL,(char*)&end,
	 (char*)"Last element to be considered (<0: go to the end)"},
        {(char*)"-iter",TK_ARGV_INT,(char*)NULL,(char*)&iter,
         (char*)"Number of iterations [default is 3]"}
      };
      
      if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,TK_ARGV_NO_LEFTOVERS))!=TCL_OK) {
        return error;
      }
      
      if (!beamline_name && !beam_name) {
        Tcl_SetResult(interp,"SetReferenceEnergy: you must provide both  '-beamline' and '-beam'", TCL_VOLATILE);
        return TCL_ERROR;
      }
      
      BEAMLINE *beamline_ptr = get_beamline(beamline_name);
      if (!beamline_ptr) beamline_ptr = inter_data.beamline;

      BEAM *beam_ptr;
      if (beam_name) {
        beam_ptr=get_beam(beam_name);
      } else {
        beam_ptr=inter_data.bunch;
      }
      
      void (*survey)(BEAMLINE*);
      if ((survey=(void (*)(BEAMLINE *))survey_find(survey_name))==NULL) {
        check_error(interp,argv[0],error);
        Tcl_AppendResult(interp,"-survey ",survey_name," not found\n",NULL);
      }

      start = std::max(0, std::min(start, beamline_ptr->n_elements));
      
      if (end<0)
	end = beamline_ptr->n_elements;
      else
	end = std::min(end, beamline_ptr->n_elements);
      
      BEAM *tb=bunch_remake(beam_ptr);
      survey(beamline_ptr);
      for (int i=0; i<iter; i++) {
	beam_copy(beam_ptr,tb);
	ELEMENT **element=beamline_ptr->element;
	inter_data.bunch=tb;
	double beam_energy0 = get_beam_energy(tb);
	for (;start<end;start++) {
	  element[start]->track_0(tb);
	  double beam_energy = get_beam_energy(tb);
	  double old_reference_energy = element[start]->get_ref_energy();
	  double new_reference_energy = (beam_energy0 + beam_energy) / 2;
	  if (QUADRUPOLE *quad_ptr = element[start]->quad_ptr()) {
	    double strength = quad_ptr->get_strength();
	    quad_ptr->set_strength(strength / old_reference_energy * new_reference_energy);
	  } else if (QUADBPM *quadbpm_ptr = element[start]->quadbpm_ptr()) {
	    double strength = quadbpm_ptr->get_strength();
	    quadbpm_ptr->set_strength(strength / old_reference_energy * new_reference_energy);
	  } else if (MULTIPOLE *multi_ptr = element[start]->multipole_ptr()) {
	    std::complex<double> strength = multi_ptr->get_strength();
	    multi_ptr->set_strength(strength / old_reference_energy * new_reference_energy);
	  } else if (SBEND *sbend_ptr = element[start]->sbend_ptr()) {
	    sbend_ptr->set_ref_energy(new_reference_energy);
	  } 
	  element[start]->set_ref_energy(new_reference_energy);  
	  beam_energy0 = beam_energy;
	}
      }
      beam_delete(tb);
    }
    return TCL_OK;
  }

  /***********************************************************************
   **                        NameNumberList                              **
   ***********************************************************************/

  int tk_NameNumberList(ClientData /*clientdata*/,Tcl_Interp *interp,int argc,
  			char *argv[])
  {
    char *beamlinename=NULL;
    Tk_ArgvInfo table[] = {
      {(char*)"-beamline",TK_ARGV_STRING,(char*)NULL,(char*)&beamlinename,
       (char*)"Name of the beamline to be used"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK) {
      return error;
    }
    BEAMLINE *beamline = get_beamline(beamlinename);
    if (!beamline) beamline = inter_data.beamline;
    ELEMENT **element=beamline->element;
    char buffer[100];
    for (int i=0;i<beamline->n_elements;i++) {
      for (int arg=1;arg<argc;arg++) {
        if (fnmatch(argv[arg],element[i]->get_name(),0)==0) {
          snprintf(buffer,100,"%d",i);
          Tcl_AppendElement(interp,buffer);
          break;	  
        }
      }
    }
    return TCL_OK;
  }
  
  /**********************************************************************/
  /*                          ElementGetName                          */
  /**********************************************************************/
  
  int tk_ElementGetName(ClientData /*clientdata*/,Tcl_Interp *interp,
			int argc,char *argv[])
  {
    int error,j=-1;
    Tk_ArgvInfo table[]={
      {(char*)NULL,TK_ARGV_HELP,(char*)NULL,(char*)NULL,
       (char*)"This command returns the name of the element number specified."},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    if (argc!=2){
      Tcl_AppendResult(interp,"Error in \n",argv[0],"\n",
		       "Usage: ",argv[0]," number\n",NULL);
      return TCL_ERROR;
    }
    if (error=Tcl_GetInt(interp,argv[1],&j)) return error;
    if (!inter_data.beamline->element_in_range(j)){return TCL_ERROR;}
    char buf[200];
    snprintf(buf,200,inter_data.beamline->element[j]->get_name());
    Tcl_SetResult(interp,buf,TCL_VOLATILE);
    return TCL_OK;  
  }


  /**********************************************************************/
  /*                          ElementSetAttributes                      */
  /**********************************************************************/
  
  int tk_ElementSetAttributes(ClientData /*clientdata*/,Tcl_Interp *interp,
			      int argc,char *argv[])
  {
    Tk_ArgvInfo table[]={
      {(char*)NULL,TK_ARGV_HELP,(char*)NULL,(char*)NULL,
       (char*)"Set one or more attributes of an element (or a list of elements)."},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    if (argc<2){
      Tcl_AppendResult(interp,"Error in \n",argv[0],"\n",
		       "Usage: ",argv[0]," { 1 2 3 } -length 123.0 -synrad 1\n",NULL);
      return TCL_ERROR;
    }
    for (int j=2;j<argc;j+=2) {
      if (argv[j][0]!='-') {
        Tcl_AppendResult(interp,"Error in \n",argv[0],"\n", "invalid attribute name '", argv[j], "'.\n",NULL);
        return TCL_ERROR;
      }
    }
    int nelem;
    char **velem;
    if (error=Tcl_SplitList(interp,argv[1],&nelem,&velem)) return error;
    char *_argv[argc-1];
    _argv[0] = NULL;
    for (int i=0;i<nelem;i++) {
      int elem;
      Tcl_GetInt(interp, velem[i], &elem);    
      if (!inter_data.beamline->element_in_range(elem)){return TCL_ERROR;}
      for (int j=1;j<argc-1;j++)
        _argv[j] = strdup(argv[j+1]);
      int _argc = argc-1;
      inter_data.beamline->element[elem]->set_attributes(_argc,_argv);
      for (int j=1;j<argc-1;j++)
        free(_argv[j]);
    }
    return TCL_OK;  
  }

  /**********************************************************************/
  /*                          ElementGetAttribute                       */
  /**********************************************************************/

  int tk_ElementGetAttribute(ClientData /*clientdata*/,Tcl_Interp *interp,
			     int argc,char *argv[])
  {
    Tk_ArgvInfo table[]={
      {(char*)NULL,TK_ARGV_HELP,(char*)NULL,(char*)NULL,
       (char*)"Get an attribute of an element (or a list of elements)."},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    if (argc<2){
      Tcl_AppendResult(interp,"Error in \n",argv[0],"\n",
		       "Usage: ",argv[0]," { 1 2 3 } -length\n",NULL);
      return TCL_ERROR;
    }
    if (argv[2][0]!='-') {
      Tcl_AppendResult(interp,"Error in ",argv[0],": invalid attribute name '", argv[2], "'. Attribute name should start with '-' (e.g. '-length')\n",NULL);
      return TCL_ERROR;
    }
    int nelem;
    char **velem;
    if (error=Tcl_SplitList(interp,argv[1],&nelem,&velem)) return error;
    for (int i=0;i<nelem;i++) {
      int elem;
      Tcl_GetInt(interp, velem[i], &elem);    
      if (!inter_data.beamline->element_in_range(elem)){return TCL_ERROR;}
      std::string str = inter_data.beamline->element[elem]->get_attribute_as_string(argv[2]+1);
      Tcl_SetResult(interp,const_cast<char*>(str.c_str()),TCL_VOLATILE);
    }
    return TCL_OK;  
  }
  
  /**********************************************************************/
  /*                          VaryCorrector                             */
  /**********************************************************************/
  
  int tk_VaryCorrector(ClientData /*clientdata*/,Tcl_Interp *interp,
		       int argc,char *argv[])
  {
    double x=0.0,y=0.0;
    Tk_ArgvInfo table[]={
      {(char*)"-x",TK_ARGV_FLOAT,(char*)NULL,
       (char*)&x,
       (char*)"horizontal shift"},
      {(char*)"-y",TK_ARGV_FLOAT,(char*)NULL,
       (char*)&y,
       (char*)"vertical shift"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    if (argc<2){
      Tcl_SetResult(interp,"Not enough arguments to <VaryCorrector>",
		    TCL_VOLATILE);
      return TCL_ERROR;
    }
    int error,n;
    if (error=Tcl_GetInt(interp,argv[1],&n)) return error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    if (argc!=2){
      Tcl_SetResult(interp,"Too many arguments to <VaryCorrector>",
		    TCL_VOLATILE);
      return TCL_ERROR;
    }
    if (x!=0.0)	inter_data.beamline->element[n]->correct_x(x);
    if (y!=0.0)	inter_data.beamline->element[n]->correct_y(y);
    return TCL_OK;
  }

  /**********************************************************************/
  /*                       ElementEnableHCorrector                      */
  /**********************************************************************/
  
  int tk_ElementEnableHCorrector(ClientData /*clientdata*/,Tcl_Interp *interp, int argc,char *argv[])
  {
    char *leverage=NULL;
    double step_size=0.0;
    Tk_ArgvInfo table[]={
      {(char*)"-leverage",TK_ARGV_STRING,(char*)NULL,
       (char*)&leverage,
       (char*)"Correction leverage [STRING]"},
      {(char*)"-step_size",TK_ARGV_FLOAT,(char*)NULL,
       (char*)&step_size,
       (char*)"Corrector step size"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    if (argc<2){
      Tcl_SetResult(interp,"Not enough arguments to <ElementEnableHCorrector>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    int error,n;
    if (error=Tcl_GetInt(interp,argv[1],&n)) return error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    if (argc!=2){
      Tcl_SetResult(interp,"Too many arguments to <ElementEnableHCorrector>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    if (leverage==NULL) {
      Tcl_SetResult(interp,"You must choose a leverage", TCL_VOLATILE);
      return TCL_ERROR;
    }
    inter_data.beamline->element[n]->enable_hcorr(leverage, step_size);
    return TCL_OK;
  }

  /**********************************************************************/
  /*                       ElementEnableVCorrector                      */
  /**********************************************************************/
  
  int tk_ElementEnableVCorrector(ClientData /*clientdata*/,Tcl_Interp *interp, int argc,char *argv[])
  {
    char *leverage=NULL;
    double step_size=0.0;
    Tk_ArgvInfo table[]={
      {(char*)"-leverage",TK_ARGV_STRING,(char*)NULL,
       (char*)&leverage,
       (char*)"Correction leverage [STRING]"},
      {(char*)"-step_size",TK_ARGV_FLOAT,(char*)NULL,
       (char*)&step_size,
       (char*)"Corrector step size"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    if (argc<2){
      Tcl_SetResult(interp,"Not enough arguments to <ElementEnableVCorrector>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    int error,n;
    if ((error=Tcl_GetInt(interp,argv[1],&n))) return error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    if (argc!=2){
      Tcl_SetResult(interp,"Too many arguments to <ElementEnableVCorrector>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    if (leverage==NULL) {
      Tcl_SetResult(interp,"You must choose a leverage", TCL_VOLATILE);
      return TCL_ERROR;
    }
    inter_data.beamline->element[n]->enable_vcorr(leverage, step_size);
    return TCL_OK;
  }

  /**********************************************************************/
  /*                       ElementDisableHCorrector                     */
  /**********************************************************************/
  
  int tk_ElementDisableHCorrector(ClientData /*clientdata*/,Tcl_Interp *interp, int argc,char *argv[])
  {
    if (argc<2){
      Tcl_SetResult(interp,"Not enough arguments to <ElementDisableHCorrector>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    int error,n;
    if ((error=Tcl_GetInt(interp,argv[1],&n))) return error;
    inter_data.beamline->element[n]->disable_hcorr();
    return TCL_OK;
  }

  /**********************************************************************/
  /*                       ElementDisableVCorrector                     */
  /**********************************************************************/
  
  int tk_ElementDisableVCorrector(ClientData /*clientdata*/,Tcl_Interp *interp, int argc,char *argv[])
  {
    if (argc<2){
      Tcl_SetResult(interp,"Not enough arguments to <ElementDisableVCorrector>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    int error,n;
    if ((error=Tcl_GetInt(interp,argv[1],&n))) return error;
    inter_data.beamline->element[n]->disable_vcorr();
    return TCL_OK;
  }

  int tk_DynamicElement(ClientData /*clientdata*/,Tcl_Interp *interp,int argc, char *argv[])
  {
    Tk_ArgvInfo table[]={
      {(char*)NULL,TK_ARGV_HELP,(char*)NULL,(char*)NULL,
       (char*)"This command places an external \"dynamically loadable\" element in the lattice."},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    try {
      DLELEMENT *dl = new DLELEMENT(argc, argv);
      inter_data.girder->add_element(dl);
    } catch (const char *err ) {
      placet_cout << ERROR << "DynamicElement error: " << err << endmsg;
      exit(1);    
    }
    return TCL_OK;
  }
  
  int tk_Verbosity(ClientData /*clientdata*/,Tcl_Interp *interp,int argc, char *argv[])
  {
    Tk_ArgvInfo table[] = {
      {(char*)NULL,TK_ARGV_HELP,(char*)NULL,(char*)NULL,
       (char*)"This command sets or returns PLACET's level of verbosity [0,6] (default: 3)."},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK) {
      return error;
    }
    if (argc==1) {
      char buffer[8];
      snprintf(buffer,8,"%d",get_verbosity());
      Tcl_AppendElement(interp,buffer);
    } else {
      set_verbosity((VERBOSITY)std::max(0,std::min(atoi(argv[1]),(int)DEBUG)));
    }      
    return TCL_OK;  
  }

  int tk_BeamlineSaveLaTeX_Picture(ClientData /*clientdata*/,Tcl_Interp *interp,int argc, char *argv[])
  {
    char *filename=NULL;
    char *beamlinename=NULL;
    char *unitlength="2mm";
    double sizex=30, sizey=20;
    Tk_ArgvInfo table[]={
      {(char*)"-beamline",TK_ARGV_STRING,(char*)NULL,(char*)&beamlinename,
       (char*)"Name of the beamline to be saved"},
      {(char*)"-file",TK_ARGV_STRING,(char*)NULL,(char*)&filename,
       (char*)"Name of the output file"},	
      {(char*)"-unitlength",TK_ARGV_STRING,(char*)NULL,(char*)&unitlength,
       (char*)"Unit length for the picture [default=\"2mm\"]"},	
      {(char*)"-sizex",TK_ARGV_FLOAT,(char*)NULL,(char*)&sizex,
       (char*)"Horizontal size of the picture [UNIT LENGTH]"},
      {(char*)"-sizey",TK_ARGV_FLOAT,(char*)NULL,(char*)&sizex,
       (char*)"Vertical size of the picture [UNIT LENGTH]"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }

    BEAMLINE *beamline = get_beamline(beamlinename);
    if (!beamline) beamline = inter_data.beamline;

    if (filename) {
      LaTeX_Picture file(filename,sizex,sizey,unitlength);
      file << *beamline;
    } else {
      LaTeX_Picture file(sizex,sizey,unitlength);
      file << *beamline;
    }
    return TCL_OK;  
  }

  int tk_BeamlineSaveFootprint(ClientData /*clientdata*/,Tcl_Interp *interp,int argc, char *argv[])
  {
    int start=0,end=-1;
    char *filename=NULL;
    char *beamlinename=NULL;
    Tk_ArgvInfo table[]={
      {(char*)"-beamline",TK_ARGV_STRING,(char*)NULL,(char*)&beamlinename,
       (char*)"Name of the beamline to be saved"},
      {(char*)"-file",TK_ARGV_STRING,(char*)NULL,(char*)&filename,
       (char*)"Name of the output file"},	
      {(char*)"-start",TK_ARGV_INT,(char*)NULL,(char*)&start,
       (char*)"Element to start from [INTEGER]"},
      {(char*)"-end",TK_ARGV_INT,(char*)NULL,(char*)&end,
       (char*)"Element to end to (included) [INTEGER]"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    BEAMLINE *beamline = get_beamline(beamlinename);
    if (!beamline) beamline = inter_data.beamline;

    if (filename) {
      std::ofstream file(filename);
      beamline->save_footprint(file, start, end);
    } else {
      beamline->save_footprint(std::cout, start, end);
    }
    return TCL_OK;  
  }

  
  int tk_LinkElements(ClientData /*clientdata*/,Tcl_Interp *interp, int argc,char *argv[])
  {
    char *beamlinename=NULL;
    Tk_ArgvInfo table[]={
      {(char*)"-beamline",TK_ARGV_STRING,(char*)NULL,(char*)&beamlinename,
       (char*)"Name of the beamline to be saved"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    int error;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK){
      return error;
    }
    BEAMLINE *beamline = get_beamline(beamlinename);
    if (!beamline) beamline = inter_data.beamline;
    if (argc<3){
      Tcl_SetResult(interp, "Not enough arguments to <LinkElements>", TCL_VOLATILE);
      return TCL_ERROR;
    }
    int n, link_to_n;
    if ((error=Tcl_GetInt(interp,argv[1],&n))) return error;
    if ((error=Tcl_GetInt(interp,argv[2],&link_to_n))) return error;
    LINK *link = new LINK(beamline->element[n]);
    link->next = beamline->element[link_to_n]->next;
    GIRDER *girder=beamline->first;
    while (girder!=NULL) {
      ELEMENT *element=girder->element();
      while (element!=NULL) {
        if (element->next==beamline->element[link_to_n]) {
          element->next=link;
        }
        element=element->next;
      }
      if (girder->element()==beamline->element[link_to_n]) {
        girder->set_first_element(link);
      }
      girder=girder->next();
    }
    delete beamline->element[link_to_n];
    beamline->unset();
    // set beamline, output suppressed
    beamline->set(NULL);
    return TCL_OK;
  }
  
  int tk_GetTransferMatrix(ClientData /*clientdata*/,Tcl_Interp *interp,int argc, char *argv[])
  {
    int error;
    int start=0,end=-1;
    char *beamlinename=NULL;
    Tk_ArgvInfo table[]={ 
      {(char*)"-beamline",TK_ARGV_STRING,(char*)NULL,(char*)&beamlinename,
       (char*)"Name of the beamline to be saved"},
      {(char*)"-start",TK_ARGV_INT,(char*)NULL,(char*)&start,
       (char*)"Element to start from [INTEGER]"},
      {(char*)"-end",TK_ARGV_INT,(char*)NULL,(char*)&end,
       (char*)"Element to end to [INTEGER]"},
      {(char*)NULL,TK_ARGV_HELP,(char*)NULL,(char*)NULL,
       (char*)"This command returns the transfer matrix from element 'start' to element 'end'"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,0))!=TCL_OK) {
      return error;
    }
    BEAMLINE *beamline = get_beamline(beamlinename);
    if (!beamline) beamline = inter_data.beamline;
    if (end==-1 || end>=beamline->n_elements)
      end=beamline->n_elements-1;
    if (start > end) {
      Tcl_SetResult(interp,"Parameter 'start' must be smaller or equal to parameter 'end'", TCL_VOLATILE);
      return TCL_ERROR;
    }
    Matrix<6,6> R = Identity<6,6>();
    for (int i=start; i<=end; i++) {
      ELEMENT *element = beamline->element[i];
      R *= element->get_transfer_matrix_6d();
    }
    Tcl_ResetResult(interp);
    for (int i=0;i<6;i++) {
      char line[1024]={'\0'};
      for (int j=0;j<6;j++) {
        char val[128];
	if (j==0)
	  snprintf(val, 128, "%g", R[i][j]);
	else
	  snprintf(val, 128, " %g", R[i][j]);
	strncat(line, val, 128); // safe since 6*128 < 1024
      }
      Tcl_AppendElement(interp, line);  
    }
    return TCL_OK;
  }

  int tk_Version(ClientData /*clientdata*/,Tcl_Interp *interp,int argc, char *argv[])
  {
#ifdef PLACET_SVN_REVISION
    Tcl_SetResult(interp, "PLACET Version No " PACKAGE_VERSION " (SVN revision r" PLACET_SVN_REVISION ")", TCL_VOLATILE);
#else
    Tcl_SetResult(interp, "PLACET Version No " PACKAGE_VERSION, TCL_VOLATILE);
#endif
    return TCL_OK;
  }

#if MPI_MODULE == 1
  int tk_MPI_TestNoCorrection(ClientData /*clientdata*/, Tcl_Interp *interp, int argc, char *argv[])
  {
    int error;
    char *file_name=NULL;
    char *beamname=NULL,*survey_name,sv[]="Zero";
    char *format = default_format;
    char *mpirun_cmd="mpirun -np 4";
    void (*survey)(BEAMLINE*);
    double bpm_res=0.0;
    BEAM *beam=NULL;
    int iter=1;
    Tk_ArgvInfo table[]={
      {(char*)"-machines",TK_ARGV_INT,(char*)NULL,(char*)&iter,
       (char*)"Number of machines to simulate"},
      {(char*)"-beam",TK_ARGV_STRING,(char*)NULL,(char*)&beamname,
       (char*)"Name of the beam to be used for tracking"},
      {(char*)"-survey",TK_ARGV_STRING,(char*)NULL,(char*)&survey_name,
       (char*)"Type of prealignment survey to be used"},
      {(char*)"-emitt_file",TK_ARGV_STRING,(char*)NULL,(char*)&file_name,
       (char*)"Filename for the results defaults to NULL (no output)"},
      {(char*)"-format",TK_ARGV_STRING,(char*)NULL,(char*)&format,
       (char*)"emittance file format"},
      {(char*)"-bpm_res",TK_ARGV_FLOAT,(char*)NULL,(char*)&bpm_res,
       (char*)"BPM resolution"},
      {(char*)"-mpirun",TK_ARGV_STRING,(char*)NULL,(char*)&mpirun_cmd,
       (char*)"Command for launching MPI"},
      {(char*)NULL,TK_ARGV_END,(char*)NULL,(char*)NULL,(char*)NULL}
    };

    survey_name=sv;
    if ((error=Tk_ParseArgv(interp,NULL,&argc,argv,table,TK_ARGV_NO_LEFTOVERS))!=TCL_OK){
      return error;
    }

    beamline_survey_hook_interp=interp;
    if (check_beamline_ready(interp,argv[0],error)) return TCL_ERROR;

    if (argc!=1){
      Tcl_SetResult(interp,"Too many arguments to TestNoCorrection",TCL_VOLATILE);
      return TCL_ERROR;
    }
    if (beamname!=NULL){
      beam=get_beam(beamname);
    } else {
      check_error(interp,argv[0],error);
      Tcl_AppendResult(interp,"-beam must be defined\n",NULL);
      error=TCL_ERROR;
    }
    if (!(survey=(void (*)(BEAMLINE *))survey_find(survey_name))) {
      check_error(interp,argv[0],error);
      Tcl_AppendResult(interp,"-survey ",survey_name," not found\n",NULL);
    }
    if (error)
      return error;

    test_no_correction_parallel(inter_data.beamline, beam, iter, survey, file_name, format, mpirun_cmd);
    return TCL_OK;
  }
#endif
}

namespace Andrea {
  
  bool UseDispersionFreeEmittance = false;
  
  double ReferenceEnergy = -1;
    
#define sigma_yy sigma
  EMITT_MINIMIZATION(y)
#ifdef TWODIM
  EMITT_MINIMIZATION(x)
#endif
#undef sigma_yy

    Tcl_Interp *Placet_interp;

  int Init(Tcl_Interp *interp )
  {
    Placet_CreateCommand(interp, "DispersionFreeEmittance", &tk_DispersionFreeEmittance, NULL, NULL);
    Placet_CreateCommand(interp, "SetReferenceEnergy", &tk_SetReferenceEnergy, NULL, NULL);
    Placet_CreateCommand(interp, "NameNumberList", &tk_NameNumberList, NULL, NULL);
    Placet_CreateCommand(interp, "ElementGetName", &tk_ElementGetName, NULL, NULL);
    Placet_CreateCommand(interp, "VaryCorrector", tk_VaryCorrector, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "ElementEnableHCorrector", tk_ElementEnableHCorrector, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "ElementEnableVCorrector", tk_ElementEnableVCorrector, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "ElementDisableHCorrector", tk_ElementDisableHCorrector, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "ElementDisableVCorrector", tk_ElementDisableVCorrector, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "ElementSetAttributes", tk_ElementSetAttributes, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "ElementGetAttribute", tk_ElementGetAttribute, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "ParticlesToSlices", &tk_ParticlesToSlices, NULL, NULL);
    
    Placet_CreateCommand(interp, "KalmanCorrectionInit", Tcl_KalmanCorrectionInit, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    
    Placet_CreateCommand(interp, "CollimatorNumberList", tk_CollimatorNumberList, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

#if OCTAVE_INTERFACE == 1
    Placet_CreateCommand(interp, "octave", tk_Octave, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "Octave", tk_Octave, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#endif

#if PYTHON_INTERFACE == 1
    Placet_CreateCommand(interp, "python", tk_Python, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "Python", tk_Python, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#endif

    Placet_CreateCommand(interp, "DynamicElement", tk_DynamicElement, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "CrabCavity", tk_CrabCavity, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "Verbosity", tk_Verbosity, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "BeamlineSaveAML", tk_BeamlineSaveAML, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "BeamlineSaveFootprint", tk_BeamlineSaveFootprint, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "BeamlineSaveLaTeX_Picture", tk_BeamlineSaveLaTeX_Picture, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "LinkElements", tk_LinkElements, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "GetTransferMatrix", tk_GetTransferMatrix, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#ifdef HTGEN
    Placet_CreateCommand(interp, "BeamReadHalo",tk_BeamReadHalo, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "BeamWriteHalo",tk_BeamWriteHalo, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "BeamClearHalo",tk_BeamClearHalo, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#endif
    Placet_CreateCommand(interp, "Manual", tk_Manual, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

    Placet_CreateCommand(interp, "version", tk_Version, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
    Placet_CreateCommand(interp, "Version", tk_Version, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);

#if MPI_MODULE == 1
    Placet_CreateCommand(interp, "MPI_Begin", &tk_BeginParallel, NULL, NULL);
    Placet_CreateCommand(interp, "MPI_End", &tk_EndParallel, NULL, NULL);
    Placet_CreateCommand(interp, "MPI_TestNoCorrection", tk_MPI_TestNoCorrection, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
#endif

    return 0;
  }
  
}
