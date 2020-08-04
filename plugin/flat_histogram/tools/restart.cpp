#include <fstream>
#include "feasst.h"

feasst::ArgumentParse parser(
"Restart simulations from checkpoint files.\n\n"

"options:\n"
"-h : print this documentation to screen.\n"
"-s : include this option for serial simulations.\n"
"-f : checkpoint file name to append (default: clone).\n"
"-p : prepend this string onto the file name (default: .fst)\n"
"-n : number of clones (default: 12)\n"
"-m : minimum index (default: 0)\n\n"

"For example, to restart a parallel simulation with the following files:\n"
"(clone0.fst, clone1.fst, clone2.fst),\n\n"

"./restart -f clone -p .fst -n 3\n");

int main(int argc, char ** argv) {
  parser.parse(argc, argv);
  std::string file_name = parser.get("-f", "clone");

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
    parser.get_int("-n", "12"),
    parser.get_int("-m", "0"),
    parser.get("-p", ".fst"));
  clones->initialize_and_run_until_complete();
}

