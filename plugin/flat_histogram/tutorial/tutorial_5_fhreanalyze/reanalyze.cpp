#include "feasst.h"

static feasst::ArgumentParse args("reanalyze dump files.", {
  {"--first_proc", "id of the first processor", "3"},
  {"--last_proc", "id of the last processor", "3"},
});

int main(int argc, char ** argv) {
  std::cout << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;
  for (int proc = args.get_int("--first_proc");
       proc <= args.get_int("--last_proc");
       ++proc) {
    feasst::MonteCarlo mc;
    feasst::Checkpoint({{"file_name", "../checkpoint" + feasst::str(proc) + ".fst"}}).read(&mc);
    feasst::FlatHistogram fh(mc.criteria());
    const feasst::Histogram& hist = fh.macrostate().histogram();
    auto tm = feasst::MakeTransitionMatrix({{"min_sweeps", "1"},
                                            {"num_blocks", "0"}});
    tm->resize(hist);
    tm->read_dump_file("../dump" + feasst::str(proc) + ".txt");
    tm->infrequent_update();
    for (const double val : tm->ln_prob().values()) {
      std::cout << val << std::endl;
    }
  }
}
