#include <stdio.h>
#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
#include "feasst/include/feasst.h"

using namespace feasst;

// HWH no need for parse_args to be a function
void parse_args(const std::string major, std::stringstream& ss, argtype * args, argtype * variables, bool * assign_to_list) {
}

std::pair<std::string, argtype> parse_line(const std::string line,
  argtype * variables,
  bool * assign_to_list) {
  std::stringstream ss(line);
  std::string major;
  ss >> major;
  argtype args;
  while(!ss.eof()) {
    std::string minor, value;
    ss >> minor >> value;
    ASSERT(!value.empty(), "Error parsing text file on line: \"" << ss.str()
      << "\". Line syntax typically requires an odd number of space-separated "
      << "strings (e.g., Object key0 value0 ... keyN valueN."
      << " This error typically occurs when one of a key/value pair is missing.");
    DEBUG("major " << major << " minor " << minor << " value " << value);
    if (major == "set_variable") {
      DEBUG("setting variable");
      (*variables)[minor] = value;
      *assign_to_list = false;
    } else if (variables->count(value) > 0) {
      DEBUG("using variable");
      args[minor] = (*variables)[value];
    } else {
      DEBUG("no variable: " << value << " sz " << value.size());
      if (value.size() > 7) {
        DEBUG(value.substr(0, 7));
        if (value.substr(0, 7) == "/feasst") {
          DEBUG("replaced: " << value);
          value.replace(0, 7, install_dir());
        }
      }
      args[minor] = value;
    }
  }
  return std::pair<std::string, argtype>(major, args);
}

arglist parse_mc(argtype variables = argtype()) {
  std::string line;
  arglist list;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line[0] != '#') {
      bool assign_to_list = true;
      std::pair<std::string, argtype> line_pair = parse_line(line, &variables, &assign_to_list);
      if (assign_to_list) {
        list.push_back(line_pair);
      }
    }
  }
  return list;
}

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
  //INFO("cm args: " << str(line_pair.second));
  CollectionMatrixSplice cm(line_pair.second);

  // read next line for Window
  std::getline(std::cin, line);
  ASSERT(line.substr(0, 6) == "Window",
    "CollectionMatrixSplice must be followed with Window. Instead, " << line);
  line_pair = parse_line(line, &variables, &assign_to_list);
  std::shared_ptr<Window> window;
  if (line_pair.first == "WindowExponential") {
    window = std::make_shared<WindowExponential>(line_pair.second);
  } else {
    FATAL("unrecognized: " << line_pair.first);
  }
  ASSERT(window->overlap() == 0, "CollectionMatrixSplice requires 0 overlap");

  // read next line for Checkpoint
  std::getline(std::cin, line);
  ASSERT(line.substr(0, 10) == "Checkpoint",
    "CollectionMatrixSplice must be followed with Window then Checkpoint. Instead, " << line);
  line_pair = parse_line(line, &variables, &assign_to_list);
  cm.set(MakeCheckpoint(line_pair.second));

  arglist list = parse_mc();
  std::vector<int> complete(window->num(), 0);
  #ifdef _OPENMP
  #pragma omp parallel
  {
    const int thread = omp_get_thread_num();
    ASSERT(thread < window->num(), "thread: " << thread << " should be less "
      << "than the number of windows: " << window->num());
  #else // _OPENMP
  for (int thread = 0; thread < window->num(); ++thread) {
  #endif // _OPENMP
    arglist list2 = list;
    replace_value("[soft_macro_min]", str(window->boundaries()[thread][0]), &list2);
    replace_value("[soft_macro_max]", str(window->boundaries()[thread][1]), &list2);
    replace_in_value("[sim_index]", str(thread), &list2);
    auto mc = std::make_shared<MonteCarlo>(list2);
    #ifdef _OPENMP
    #pragma omp critical
    #endif // _OPENMP
    {
      cm.set(thread, mc);
    }
    complete[thread] = 1;
    while (std::accumulate(complete.begin(), complete.end(), 0) < window->num()) {
      cm.get_clone(thread)->attempt(10);
    }
  }
  cm.run_until_all_are_complete();
}

/**
  Usage: ./fst < file.txt

  Syntax: MAJOR MINOR_0 VALUE_0 ... MINOR_N VALUE_N

  "set_variable name value" will replace any use of name in subsequent values.

  Any value beginning with "/feasst" will have that beginning replaced with
  feasst::install_dir().
 */
int main() {
  std::cout << "# FEASST version: " << version() << std::endl;
  std::string line;
  std::getline(std::cin, line);

  // Find the first non-empty and non-comment line in the script
  int num_lines = 0;
  while (line.empty() || line[0] == '#') {
    std::getline(std::cin, line);
    ++num_lines;
    ASSERT(num_lines < 1e4, "Improperly formated input");
  }

  if (line == "MonteCarlo") {
    arglist list = parse_mc();
    auto mc = std::make_shared<MonteCarlo>(list);
  } else if (line.substr(0, 22) == "CollectionMatrixSplice") {
    parse_cm(line);
  } else {
    FATAL("As currently implemented, all FEASST input text files must begin "
      << "with \"MonteCarlo\" or \"CollectionMatrixSplice\", but the first "
      << "readable line in this file is: " << line);
  }
  return 0;
}

