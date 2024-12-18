#include <memory>
#include <iostream>
#include <string.h>
#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
#include "feasst/include/feasst.h"

using namespace feasst;

// Parse the line containing CollectionMatrixSplice and the following line
// beginning "Window"
// consider moving this to CollectionMatrixSplice
// - read arglist like in MC, but including the first line containing CollectionMatrixSplice
// - parse arglist through CollectionMatrixSplice constructor
void parse_cm(std::string line) {
  // parse and obtain Window arguments here, as well as some additional CollectionMatrixSplice arguments.
  // use the windows to set custom variables here for soft_min and soft_max.
  // for each proc, input to bounds as variables into parse_mc to initialize clones argument
  // finally, run CollectionMatrixSplice.
  // restart should try { MonteCarlo } and, if fail, try CollectionMatrixSplice, wish some kind of resume function.
  argtype variables;
  bool assign_to_list;
  std::pair<std::string, argtype> line_pair = parse_line(line, &variables, &assign_to_list);
  //DEBUG("cm args: " << str(line_pair.second));
  CollectionMatrixSplice cm(line_pair.second);

  // read next line for Window
  std::getline(std::cin, line);
  ASSERT(line.substr(0, 6) == "Window",
    "CollectionMatrixSplice must be followed with Window. Instead, " << line);
  line_pair = parse_line(line, &variables, &assign_to_list);
  std::shared_ptr<Window> window;
  if (line_pair.first == "WindowExponential") {
    window = std::make_shared<WindowExponential>(line_pair.second);
  } else if (line_pair.first == "WindowCustom") {
    window = std::make_shared<WindowCustom>(line_pair.second);
  } else {
    FATAL("unrecognized: " << line_pair.first);
  }
  ASSERT(window->overlap() == 0, "CollectionMatrixSplice requires 0 overlap");

  // read next line for Checkpoint
  std::getline(std::cin, line);
  ASSERT(line.substr(0, 10) == "Checkpoint",
    "CollectionMatrixSplice must be followed with Window then Checkpoint. " <<
    "Instead, the following line was given: " << line);
  line_pair = parse_line(line, &variables, &assign_to_list);
  cm.set(MakeCheckpoint(line_pair.second));

  std::vector<arglist> list = parse_mcs(std::cin);
  ASSERT(static_cast<int>(list.size()) == 1, "CollectionMatrixSplice should "
    << "have no lines beginning as \"MonteCarlo\"");
  std::vector<int> complete(window->num(), 0);
  cm.set_size(window->num());
  #ifdef _OPENMP
  #pragma omp parallel
  {
    const int num_threads = MakeThreadOMP()->num();
    ASSERT(num_threads >= window->num(),
      "asked for " << window->num() << " windows but there are only "
      << num_threads << " OMP threads");
    const int thread = omp_get_thread_num();
    if (thread < window->num()) {
      arglist list2 = list[0];
      replace_value("[soft_macro_min]", str(window->boundaries()[thread][0]), &list2);
      replace_value("[soft_macro_max]", str(window->boundaries()[thread][1]), &list2);
      replace_in_value("[sim_index]", sized_int_to_str(thread, num_threads), &list2);
      auto mc = std::make_shared<MonteCarlo>(list2);
      cm.set(thread, mc);
      complete[thread] = 1;
      while (std::accumulate(complete.begin(), complete.end(), 0) < window->num()) {
        cm.get_clone(thread)->attempt(10);
      }
    }
  }
  #else // _OPENMP
    FATAL("CollectionMatrixSplice requires OMP");
  #endif // _OPENMP
  cm.run_until_all_are_complete();
}

// Parse the line containing Prefetch
void parse_prefetch(std::string line) {
  argtype variables;
  bool assign_to_list;
  std::pair<std::string, argtype> line_pair = parse_line(line, &variables, &assign_to_list);
  std::vector<arglist> lists = parse_mcs(std::cin);
  for (auto list : lists) {
    Prefetch prefetch(line_pair.second);
    prefetch.begin(list);
  }
}

void parse_server(std::string line) {
  argtype variables;
  bool assign_to_list;
  std::pair<std::string, argtype> line_pair = parse_line(line, &variables, &assign_to_list);
  ASSERT(line_pair.first == "Server", "error");
  auto server = std::make_shared<Server>(line_pair.second);
  server->bind_listen_accept();
  bool finished = false;
  std::string sim_type;
  std::shared_ptr<MonteCarlo> mc;
  while (!finished) {
    const int size = server->receive();
    ASSERT(size >= 0, "error");
    if (strcmp(server->buffer(), "EndListen") == 0) {
      DEBUG(server->buffer());
      finished = true;
    } else {
      if (size > 0) {
        std::string line(server->buffer());
        std::pair<std::string, argtype> marg = parse_line(line, NULL, NULL);
        DEBUG(marg.first);
        DEBUG(str(marg.second));
        if (marg.first == "MonteCarlo") {
          // std::cout << "MonteCarlo" << std::endl;
          sim_type = marg.first;
          mc = std::make_shared<MonteCarlo>();
        } else {
          arglist list(1);
          list[0] = marg;
          if (sim_type.empty()) {
            FATAL(marg.first << " not yet implemented in server-client mode.");
          } else if (sim_type == "MonteCarlo") {
            mc->parse_args(&list, true); // silence std::cout
          }
          ASSERT(list.size() == 0, "Unrecognized argument: " << list.begin()->first);
        }
      }
      server->send("received");
    }
    ASSERT(size >= 0, "error");
  }
}

std::unique_ptr<MonteCarlo> restart(const std::string& filename, bool resume = true) {
  std::unique_ptr<MonteCarlo> mc;
  CollectionMatrixSplice cms;
  bool is_mc = false;
  try {
    MakeCheckpoint({{"checkpoint_file", filename}})->read(&cms);
  } catch (const CustomException& e) {
    MakeCheckpoint({{"checkpoint_file", filename}})->read_unique(mc);
    is_mc = true;
  }
  if (is_mc) {
    if (resume) {
      mc->resume();
    } else {
      mc->clear_arguments();
    }
  } else {
    cms.run_until_all_are_complete();
  }
  return mc;
}

void parse_restart(std::string line) {
  std::stringstream ss(line);
  std::string checkpoint_file;
  ss >> checkpoint_file; // reads "Restart"
  ss >> checkpoint_file;
  std::string extra;
  ss >> extra;
  DEBUG("extra " << extra);
  bool resume = true;
  if (extra == "clear_previous_arguments") {
    resume = false;
  }
  DEBUG("resume " << resume);
  std::cout << "# Restarting from file: " << checkpoint_file << std::endl;
  std::unique_ptr<MonteCarlo> mc = restart(checkpoint_file, resume);
  // after restart is complete, check for more text input
  if (mc) {
    for (arglist argl : parse_mcs(std::cin)) {
      mc->add_args(argl);
    }
    mc->resume();
  }
 // arglist list;
 // argtype variables;
//  while (std::getline(std::cin, line)) {
//    if (!line.empty() && line[0] != '#') {
//      bool assign_to_list = true;
//      std::pair<std::string, argtype> line_pair = parse_line(line, &variables, &assign_to_list);
//      if (assign_to_list) {
//        list.push_back(line_pair);
//      }
//    }
//  }
//  mc->begin(list);
}

/**
  Usage: ./fst < file.txt

  Syntax: MAJOR MINOR_0 VALUE_0 ... MINOR_N VALUE_N

  "set_variable name value" will replace any use of name in subsequent values.

  Any value beginning with "/feasst" will have that beginning replaced with
  feasst::install_dir().
 */
int main() {
  std::cout << "# Usage: ./fst < file.txt" << std::endl;
  std::cout << "FEASST version " << version() << std::endl;
  std::string line;
  std::getline(std::cin, line);

  // Find the first non-empty and non-comment line in the script
  int num_lines = 0;
  while (line.empty() || line[0] == '#') {
    std::getline(std::cin, line);
    ++num_lines;
    ASSERT(num_lines < 1e6, "Improperly formated input");
  }

  // Check for FEASST arguments and check version
  if (line.substr(0, 6) == "FEASST") {
    argtype variables;
    auto args = parse_line(line, NULL, NULL);
    const std::string read_ver = str("version", &args.second);
    if (read_ver != version()) {
      WARN("Version given in text file: " << read_ver << " is not the same " <<
           "as the executable version: " << version());
    }
    std::getline(std::cin, line);
  }

  if (line == "MonteCarlo") {
    std::cout << "MonteCarlo" << std::endl;
    std::vector<arglist> lists = parse_mcs(std::cin);
    for (auto list : lists) {
      auto mc = std::make_shared<MonteCarlo>(list);
    }
  } else if (line.substr(0, 8) == "Prefetch") {
    std::cout << "Prefetch" << std::endl;
    parse_prefetch(line);
  } else if (line.substr(0, 22) == "CollectionMatrixSplice") {
    parse_cm(line);
  } else if (line.substr(0, 6) == "Server") {
    parse_server(line);
  } else if (line.substr(0, 7) == "Restart") {
    parse_restart(line);
  } else {
    FATAL("As currently implemented, all FEASST input text files must begin "
      << "with \"MonteCarlo,\" \"CollectionMatrixSplice,\" \"Prefetch\" "
      << "or \"Server.\" The first readable line is: " << line);
  }
  return 0;
}

