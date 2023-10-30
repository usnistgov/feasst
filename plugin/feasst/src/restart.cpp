#include <stdio.h>
#include "feasst/include/feasst.h"

/**
  Usage: ./rst checkpoint.fst
 */
int main(int argc, char ** argv) {
  std::cout << "# FEASST version: " << feasst::version() << std::endl;
  ASSERT(argc == 2, "unrecognized number of arguments: " << argc);

//  // read the first word in the file.
//  std::ifstream file(file_name);
//  std::string line;
//  std::getline(file, line);
//  std::stringstream ss;
//  ss << line;
//  const std::string line

  feasst::MonteCarlo mc;
  feasst::CollectionMatrixSplice cms;
  bool is_mc = false;
  try {
    feasst::MakeCheckpoint({{"checkpoint_file", std::string(argv[1])}})->read(&cms);
  } catch (const feasst::CustomException& e) {
    feasst::MakeCheckpoint({{"checkpoint_file", std::string(argv[1])}})->read(&mc);
    is_mc = true;
  }
  if (is_mc) {
    mc.resume();
  } else {
    cms.run_until_all_are_complete();
  }
  return 0;
}
