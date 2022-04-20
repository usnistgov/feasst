#include <stdio.h>
#include "feasst/include/feasst.h"

/**
  Usage: ./rst checkpoint.fst
 */
int main(int argc, char ** argv) {
  std::cout << "# version: " << feasst::version() << std::endl;
  ASSERT(argc == 2, "unrecognized number of arguments: " << argc);
  feasst::MonteCarlo mc;
  feasst::MakeCheckpoint({{"file_name", std::string(argv[1])}})->read(&mc);
  mc.resume();
  return 0;
}
