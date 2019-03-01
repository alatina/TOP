#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <signal.h>
#include <tk.h>

static int do_wish;

static void TOP_termination_handler(int signum )
{
  exit(0);
}

static int Tcl_AppInit(Tcl_Interp *interp )
{
  /* initialisation of Tcl part */

  if (Tcl_Init(interp)){
    TOP_printf(INFO,"%s\n",Tcl_GetStringResult(interp));
  }

  if (do_wish){
    if (Tk_Init(interp)){
      TOP_printf(INFO,"%s\n",Tcl_GetStringResult(interp));
    }
  }

  TOP_Init(interp);
  return TCL_OK;
}

int main(int argc, char *argv[])
{
  // set up a termination action
  struct sigaction new_action;
  new_action.sa_handler = TOP_termination_handler;
  sigemptyset (&new_action.sa_mask);
  new_action.sa_flags = 0;
  sigaction (SIGINT, &new_action, NULL);
  // include end message at exit of TOP
  // here we go
  do_wish = 0;
  Tcl_CreateCommand(interp, "octave", tk_Octave, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp, "Octave", tk_Octave, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp, "python", tk_Python, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp, "Python", tk_Python, (ClientData)NULL,(Tcl_CmdDeleteProc*)NULL);
  if (do_wish) Tk_Main(argc,argv,&Tcl_AppInit);
  Tcl_Main(argc,argv,&Tcl_AppInit);
  exit(0);
}
