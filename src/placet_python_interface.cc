#include <vector>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#include <poll.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

//#include <readline/readline.h>
//#include <readline/history.h>

#include "placet_interface.hh"
#include "placet_api_server.hh"
#include "piped_process.hh"
#include "file_stream.hh"
#include "placet_tk.h"
#include "defaults.h"
#include "placet_prefix.hh"
#include "andrea.h"

static std::string create_script_from_args(Tcl_Interp *interp, int argc, char *argv[] );

static std::string gen_random_str(const size_t len )
{
  static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
  char str[len+1];
  for (size_t i=0; i<len; i++) {
    str[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
  }
  str[len]='\0';
  return std::string(str);
}

namespace {
  struct Python_instantiation {
    Piped_process *process_ptr;
    std::string unique_id;
    int ifiledes;
    int ofiledes;
    Python_instantiation() : process_ptr(NULL) {
      unique_id = gen_random_str(8);
      const char *argv[] = { // NULL-terminated array
	"placet-python", "-i", "-c", "import sys ; sys.ps1 = ''", NULL
      };
      process_ptr = new Piped_process(PYTHON_BIN, const_cast<char* const*>(argv));
      std::ostringstream cmd;
      cmd <<
	"import os\n"
	"sys.path.append('" << get_placet_sharepath("interfaces") << "')\n"
	"sys.path.append('" << get_placet_sharepath("python") << "')\n"
	"from numpy import *\n"
	"import placet_python as placet\n"
	"rv = placet.init(" << process_ptr->ifiledes() << "," << process_ptr->ofiledes() << ");\n";
      process_ptr->send(cmd.str().c_str());
      ifiledes = process_ptr->ifiledes();
      ofiledes = process_ptr->ofiledes();
    }
    ~Python_instantiation() {
      if (process_ptr) {
	process_ptr->send("placet.finalize()\n");
	delete process_ptr;
      }
    }
  private:
    Python_instantiation(const Python_instantiation & ) {}
  };

  std::vector<Python_instantiation*> python_instantiations;

  struct Counter {
    static int counter;
    operator int() const { return counter; }
    Counter() { ++counter; }
    ~Counter() { --counter; }
  };
  int Counter::counter = 0;
}

static void module_cleanup()
{
  for(size_t i=0; i<python_instantiations.size(); i++) {
    if (python_instantiations[i]) {
      delete python_instantiations[i];
    }
  }
}

int tk_Python(ClientData /*clientdata*/, Tcl_Interp *interp, int argc, char *argv[] )
{
  beamline_survey_hook_interp = interp;
  int return_value = TCL_OK;
  if (python_instantiations.empty()) { // first call
    atexit(module_cleanup);
  }
  Counter recursion_counter;
  python_instantiations.resize(std::max(int(recursion_counter), int(python_instantiations.size())), NULL);
  if (python_instantiations[recursion_counter-1] == NULL) // create a new instantiation
    python_instantiations[recursion_counter-1] = new Python_instantiation;
  Python_instantiation *python_ptr = python_instantiations[recursion_counter-1];
  
  bool dead = false;
  if (python_ptr&&python_ptr->ifiledes&&python_ptr->ofiledes) {
    Piped_process &python = *python_ptr->process_ptr;
    int python_stdout = python.get_stdout();
    int python_stderr = python.get_stderr();
    int ofiledes = python_ptr->ofiledes;
    int ifiledes = python_ptr->ifiledes;
    bool interactive = argc == 1;
    bool loop = true;
    if (!interactive) {
      std::string script_name = create_script_from_args(interp, argc, argv);
      if (script_name != std::string()) {
	std::ostringstream cmd;
	cmd <<
	  "execfile('" << script_name << "')\n"
	  "os.remove('" << script_name << "')\n"
	  "sys.stdout.write('" << python_ptr->unique_id << "')\n" // hook for exiting the following while(true) { ... }
	  "sys.stdout.flush()\n";
       	python.send(cmd.str().c_str());
      } else {
	loop = false;
	return_value = TCL_ERROR;
      }
    } else {
      python.send("sys.ps1 = '>>> '\n");
    }
    
    pollfd fds[4] = {
      { STDIN_FILENO,  short(interactive ? POLLIN : 0), 0 },
      { python_stdout, POLLIN, 0 },
      { python_stderr, POLLIN, 0 },
      { ofiledes,      POLLIN, 0 }
    };
    
    char buf[4096];    
    while(loop) {

      if (poll(fds, sizeof(fds) / sizeof(*fds), -1) < 0) {
	perror("poll");
	return_value = TCL_ERROR;
	break;
      }

      if (fds[0].revents) {
	int n = read(STDIN_FILENO, buf, sizeof(buf) - 1);
	if (n < 0) {
	  perror("read");
	  return_value = TCL_ERROR;
	  break;
	} else if (n==0) {
	  break;
	} else if (n > 0) {
	  if (n>=4 && strnstr(buf, "exit", 4)) {
	    break;
	  }
	  python.send(buf, n);
	}
      }

      if (fds[1].revents) {
	int n = read(python_stdout, buf, sizeof(buf) - 1);
	if (n < 0) {
	  perror("read");
	  return_value = TCL_ERROR;
	  break;
	} else if (n == 0) {
	  std::cerr << "received 0 bytes on python's stdout (python sub-process is probably not running)\n";
	  return_value = TCL_ERROR;
	  dead = true;
	  break;
	} else {
	  if (char *str = strnstr(buf, python_ptr->unique_id.c_str(), n)) {
	    write(STDOUT_FILENO, buf, str-buf);
	    break;
	  }
	  write(STDOUT_FILENO, buf, n);
	}
      }

      if (fds[2].revents) {
	int n = read(python_stderr, buf, sizeof(buf) - 1);
	if (n < 0) {
	  perror("read");
	  return_value = TCL_ERROR;
	  break;
	} else if (n == 0) {
	  std::cerr << "received 0 bytes on python's stderr (python sub-process is probably not running)\n";
	  return_value = TCL_ERROR;
	  dead = true;
	  break;
	} else {
	  if (strnstr(buf, "Error:", n)) {
	    return_value = TCL_ERROR;
	  }
	  write(STDERR_FILENO, buf, n);
	}
      }

      if (fds[3].revents) {
	int ret = placet_api_dispatcher(ofiledes, ifiledes);
	if (ret) {
	  return_value = TCL_ERROR;
	  break;
	}
      }
    }
    if (interactive && !dead) {
      python.send("sys.ps1=''\n");
    }
  }
  if (dead) {
    delete python_instantiations[recursion_counter-1];
    python_instantiations[recursion_counter-1] = NULL;
  }
  return return_value;
}

static int parse_tcl_and_write(Tcl_Interp *interp, int filedes, char *start )
{
  char *last_ptr = start;
  char *str_ptr = last_ptr;
  char last_char = '\0';
  while(*str_ptr) {
    if (*str_ptr=='$') {
      if (last_char=='\\') {
	if (str_ptr-1!=last_ptr)
	  write(filedes, last_ptr, str_ptr-last_ptr-1);
	last_ptr=str_ptr++;
      } else {
	if (str_ptr!=last_ptr)
	  write(filedes, last_ptr, str_ptr-last_ptr);
	if (const char *var = Tcl_ParseVar(interp, str_ptr, &last_ptr)) {
	  write(filedes, var, strlen(var));
	  str_ptr=last_ptr;
	} else {
	  last_ptr = NULL;
	  placet_cout << ERROR << " Tcl variable " << str_ptr << " does not exist" << endmsg;
	  return 1;
	}
      }
      last_char = *str_ptr;
    } else if (*str_ptr=='#') {
      if (last_char=='\\') {
	if (str_ptr-1!=last_ptr)
	  write(filedes, last_ptr, str_ptr-last_ptr-1);
	last_ptr=str_ptr++;
      } else {
	if (str_ptr!=last_ptr)
	  write(filedes, last_ptr, str_ptr-last_ptr);
	if ((last_ptr = strchr(str_ptr, '\n'))) {
	  str_ptr=last_ptr;
	} else {
	  return 1;
	}
      }
      last_char = *str_ptr;
    } else {
      last_char = *str_ptr++;
    }
  }
  if (str_ptr!=last_ptr)  write(filedes, last_ptr, str_ptr-last_ptr);
  return 0;
}

static std::string create_script_from_args(Tcl_Interp *interp, int argc, char *argv[] )
{
  const char *tmpdir = getenv("TMPDIR");
  if (!tmpdir) tmpdir = "/tmp/";
  std::string filename_mask = std::string(tmpdir) + std::string("/placet-python.XXXXXX");
  const int length = filename_mask.length();
  char filename[length+1];
  filename_mask.copy(filename,length,0);
  filename[length] = '\0';
  int fd = mkstemp(filename);
  if (fd != -1) {
    for (int i=1; i<argc; i++) {
      if (parse_tcl_and_write(interp, fd, argv[i])) {
	close(fd);
	remove(filename);
	return std::string();
      }
    }
    close(fd);
  } else {
    placet_cout << ERROR << ": cannot write temporary file " << std::string(filename) << endmsg;
    return std::string();
  }
  return filename;
}
