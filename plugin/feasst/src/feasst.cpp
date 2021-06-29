#include <stdio.h>
#include "feasst/include/feasst.h"

/**
  Expected usage: MAJOR MINOR_0 VALUE_0 ... MINOR_N VALUE_N

  "set_variable name value" will replace any use of name in subsequent values.

  Any value beginning with "/feasst" will have that beginning replaced with
  feasst::install_dir().
 */
int main() {
  std::cout << "version: " << feasst::version() << std::endl;
  std::string line;
  feasst::arglist list;
  feasst::argtype variables;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line[0] != '#') {
      std::stringstream ss(line);
      std::string major, minor, value;
      ss >> major;
      feasst::argtype args;
      bool assign_to_list = true;
      while(!ss.eof()) {
        ss >> minor >> value;
        DEBUG("major " << major << " minor " << minor << " value " << value);
        if (major == "set_variable") {
          DEBUG("setting variable");
          variables[minor] = value;
          assign_to_list = false;
        } else if (variables.count(value) > 0) {
          DEBUG("using variable");
          args[minor] = variables[value];
        } else {
          DEBUG("no variable: " << value << " sz " << value.size());
          if (value.size() > 7) {
            DEBUG(value.substr(0, 7));
            if (value.substr(0, 7) == "/feasst") {
              DEBUG("replaced: " << value);
              value.replace(0, 7, feasst::install_dir());
            }
          }
          args[minor] = value;
        }
      }
      if (assign_to_list) {
        list.push_back(std::pair<std::string, feasst::argtype>(major, args));
      }
    }
  }
  auto mc = std::make_shared<feasst::MonteCarlo>(list);
  return 0;
}
