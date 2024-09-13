/**
  Usage: ./rst --help
  Usage: ./rst --checkpoint_file checkpoint.fst
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "feasst/include/feasst.h"

static feasst::ArgumentParse args("Restart from checkpoint file.", {
  {"--checkpoint_file", "File name of checkpoint"},
  {"--text_file", "Optionally provide a text file with additional arguments."},
});

int main(int argc, char ** argv) {
  std::cout << "#FEASST version: " << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;
  //std::cout << "# FEASST version: " << feasst::version() << std::endl;
  //ASSERT(argc == 2, "unrecognized number of arguments: " << argc);

  ASSERT(args.option_given("--checkpoint_file"),
    "--checkpoint_file is a required argument for this executable. " <<
    "Use the --help argument for more information.");
  const std::string filename = args.get("--checkpoint_file");
  std::unique_ptr<feasst::MonteCarlo> mc;
  feasst::CollectionMatrixSplice cms;
  bool is_mc = false;
  try {
    feasst::MakeCheckpoint({{"checkpoint_file", filename}})->read(&cms);
  } catch (const feasst::CustomException& e) {
    feasst::MakeCheckpoint({{"checkpoint_file", filename}})->read_unique(mc);
    is_mc = true;
  }
  if (is_mc) {
    if (args.option_given("--text_file")) {
      const std::string filename = args.get("--text_file");
      std::cout << "# reading from: " << filename << std::endl;
      std::ifstream in(filename);
      for (feasst::arglist argl : feasst::parse_mcs(in)) {
        mc->add_args(argl);
      }
    }
    mc->resume();
  } else {
    cms.run_until_all_are_complete();
  }
  return 0;
}
