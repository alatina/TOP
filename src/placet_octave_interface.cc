#include <vector>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#include <poll.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#include "placet_interface.hh"
#include "placet_api_server.hh"
#include "piped_process.hh"
#include "file_stream.hh"
#include "placet_tk.h"
#include "defaults.h"
#include "andrea.h"
#include "placet_prefix.hh"

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
  struct Octave_instantiation {
    Piped_process *process_ptr;
    std::string unique_id;
    int ifiledes;
    int ofiledes;
    // without the use of a temporary string, argv[4] == argv[6]..
    std::string tmp_string,tmp_string2;
    Octave_instantiation() : process_ptr(NULL) {
      unique_id = gen_random_str(8);
      tmp_string =  get_placet_sharepath("interfaces");
      tmp_string2 =  get_placet_sharepath("octave");
      const char *argv[] = { // NULL-terminated array
	"placet-octave", "--quiet", "--interactive",
	  "--path", tmp_string.c_str(),
	  "--path", tmp_string2.c_str(),
	  "--eval", "_ps1=PS1('');",
	  "--persist", "--no-gui", NULL
      };
      process_ptr = new Piped_process(OCTAVE_BIN, const_cast<char* const*>(argv));
      std::ostringstream cmd;
      cmd << 
	"placet_octave;"
	"placet_octave.init(" << process_ptr->ifiledes() << ',' << process_ptr->ofiledes() << ");";
      process_ptr->send(cmd.str().c_str());
      ifiledes = process_ptr->ifiledes();
      ofiledes = process_ptr->ofiledes();
    }
    ~Octave_instantiation() {
      if (process_ptr) {
	process_ptr->send("placet_octave.finalize();\n");
	delete process_ptr;
      }
    }
  private:
    Octave_instantiation(const Octave_instantiation & ) {}
  };

  std::vector<Octave_instantiation*> octave_instantiations;

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
  for(size_t i=0; i<octave_instantiations.size(); i++) {
    if (octave_instantiations[i]) {
      delete octave_instantiations[i];
    }
  }
}

int tk_Octave(ClientData /*clientdata*/, Tcl_Interp *interp, int argc, char *argv[] )
{
  beamline_survey_hook_interp = interp;
  int return_value = TCL_OK;
  if (octave_instantiations.empty()) { // first call
    atexit(module_cleanup);
  }
  Counter recursion_counter;
  octave_instantiations.resize(std::max(int(recursion_counter), int(octave_instantiations.size())), NULL);
  if (octave_instantiations[recursion_counter-1] == NULL) // create a new instantiation
    octave_instantiations[recursion_counter-1] = new Octave_instantiation;
  Octave_instantiation *octave_ptr = octave_instantiations[recursion_counter-1];
  
  bool dead = false;
  if (octave_ptr&&octave_ptr->ifiledes&&octave_ptr->ofiledes) {
    Piped_process &octave = *octave_ptr->process_ptr;
    int octave_stdout = octave.get_stdout();
    int octave_stderr = octave.get_stderr();
    int ofiledes = octave_ptr->ofiledes;
    int ifiledes = octave_ptr->ifiledes;
    bool interactive = argc == 1;
    bool loop = true;
    if (!interactive) {
      std::string script_name = create_script_from_args(interp, argc, argv);
      if (script_name != std::string()) {
	std::ostringstream cmd;
	cmd <<
	  "source('" << script_name << "');\n"
	  "delete('" << script_name << "');\n"
	  "fputs(stdout, '" << octave_ptr->unique_id << "');\n" // hook for exiting the following while(true) { ... }
	  "fflush(stdout);\n";
       	octave.send(cmd.str().c_str());
      } else {
	loop = false;
	return_value = TCL_ERROR;
      }
    } else {
      octave.send("PS1(_ps1);\n");
    }
    
    pollfd fds[4] = {
      { STDIN_FILENO,  short(interactive ? POLLIN : 0), 0 },
      { octave_stderr, POLLIN, 0 },
      { octave_stdout, POLLIN, 0 },
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
	  octave.send(buf, n);
	}
      }

      if (fds[1].revents) {
	int n = read(octave_stderr, buf, sizeof(buf) - 1);
	if (n < 0) {
	  perror("read");
	  return_value = TCL_ERROR;
	  break;
	} else if (n == 0) {
	  std::cerr << "received 0 bytes on octave's stderr (octave sub-process is probably not running)\n";
	  return_value = TCL_ERROR;
	  dead = true;
	  break;
	} else {
	  if (strnstr(buf, "error:", n)) {
	    return_value = TCL_ERROR;
	  }
	  write(STDERR_FILENO, buf, n);
	}
      }

      if (fds[2].revents) {
	int n = read(octave_stdout, buf, sizeof(buf) - 1);
	if (n < 0) {
	  perror("read");
	  return_value = TCL_ERROR;
	  break;
	} else if (n == 0) {
	  std::cerr << "received 0 bytes on octave's stdout (octave sub-process is probably not running)\n";
	  return_value = TCL_ERROR;
	  dead = true;
	  break;
	} else {
	  if (char *str = strnstr(buf, octave_ptr->unique_id.c_str(), n)) {
	    write(STDOUT_FILENO, buf, str-buf);
	    break;
	  }
	  write(STDOUT_FILENO, buf, n);
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
      octave.send("_ps1=PS1('');\n");
    }
  }
  if (dead) {
    delete octave_instantiations[recursion_counter-1];
    octave_instantiations[recursion_counter-1] = NULL;
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
  std::string filename_mask = std::string(tmpdir) + std::string("/placet-octave.XXXXXX");
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
