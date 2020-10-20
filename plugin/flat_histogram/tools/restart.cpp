#include <fstream>
#include "feasst.h"

feasst::ArgumentParse parser("Restart simulations from checkpoint files.\n"
"For example, to restart a parallel simulation with the following files:\n"
"(clone0.fst, clone1.fst, clone2.fst),\n\n"

"./restart -f clone -p .fst -n 3", {
  {"-s", "include this option for serial simulations."},
  {"-f", "checkpoint file name to append", "clone"},
  {"-p", "prepend this string onto the file name", ".fst"},
  {"-n", "number of clones", "12"},
  {"-m", "minimum index", "0"}});

int main(int argc, char ** argv) {
  parser.parse(argc, argv);
  std::string file_name = parser.get("-f");

  // serialize option
  if (parser.option_given("-s")) {
    std::ifstream file(file_name);
    std::string line;
    std::getline(file, line);
    auto mc = feasst::MonteCarlo().deserialize(line);
    mc.run_until_complete();
    return 0;
  }

  auto clones = feasst::MakeClones(
    file_name,
    parser.get_int("-n"),
    parser.get_int("-m"),
    parser.get("-p"));
  clones->initialize_and_run_until_complete();
}

