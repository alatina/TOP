#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <signal.h>
#include <tcl.h>
#include <tk.h>

static bool do_wish = false;

int tk_Python(ClientData /*clientdata*/, Tcl_Interp *interp, int argc, char **argv )
{
  return 0;
}

int tk_Octave(ClientData /*clientdata*/, Tcl_Interp *interp, int argc, char **argv )
{
  return 0;
}

static int TOP_Init(Tcl_Interp *interp )
{
  Tcl_CreateCommand(interp, "octave", (Tcl_CmdProc *)tk_Octave, (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp, "Octave", (Tcl_CmdProc *)tk_Octave, (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp, "python", (Tcl_CmdProc *)tk_Python, (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp, "Python", (Tcl_CmdProc *)tk_Python, (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  return 0;
}

int Tcl_AppInit(Tcl_Interp *interp )
{
  /* initialisation of Tcl part */

  if (Tcl_Init(interp)){
    printf("%s\n",Tcl_GetStringResult(interp));
  }

  if (do_wish){
    if (Tk_Init(interp)){
      printf("%s\n",Tcl_GetStringResult(interp));
    }
  }

  TOP_Init(interp);
  return TCL_OK;
}

int main(int argc, char *argv[])
{
  do_wish = strncmp(argv[1], "-w", 2) == 0;
  if (do_wish) {
    --argc;
    for (int i=1; i<argc; i++) {
      argv[i] = argv[i+1];
    }
    Tk_Main(argc,argv,&Tcl_AppInit);
  }
  Tcl_Main(argc,argv,&Tcl_AppInit);
  return 0;
}
