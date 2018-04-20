/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./analyze_traj.h"

namespace feasst {

AnalyzeTRAJ::AnalyzeTRAJ(Pair *pair, const argtype &args)
  : Analyze(pair, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);

  // parse format
  format_ = argparse_.key("format").dflt("xyz").str();

  argparse_.checkAllArgsUsed();
}

AnalyzeTRAJ::AnalyzeTRAJ(
  Pair *pair,
  const char* fileName)
    : Analyze(pair, fileName) {
  defaultConstruction_();
}

void AnalyzeTRAJ::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
}

void AnalyzeTRAJ::defaultConstruction_() {
  className_.assign("AnalyzeTRAJ");
  format_.assign("xyz");
  verbose_ = 0;
}

void AnalyzeTRAJ::write() {
  if (!fileName_.empty()) {
    if (format_ == "xyz") {
      pair_->printXYZ(fileName_.c_str(), firstFlag_, space()->hash());
    #ifdef XDRFILE_H_
    } else if (format_ == "xtc") {
      stringstream ss;
      ss << fileName_ << "n" << space()->nMol();
      string mode("a");
      if (firstFlag_ == 1) {
        mode.assign("w");
      }
      pair_->printXYZ(ss.str().c_str(), 2);
      ss << ".xtc";
      XDRFILE* trjFileXDR;
      trjFileXDR = xdrfile_open(ss.str().c_str(), mode.c_str());
      space()->writeXTC(trjFileXDR);
      xdrfile_close(trjFileXDR);
    #endif  // XDRFILE_H_
    } else {
      ASSERT(0, "unrecognized format(" << format_ << ")");
    }
  }
  firstFlag_ = 0;
}

shared_ptr<AnalyzeTRAJ> makeAnalyzeTRAJ(Pair *pair, const argtype &args) {
  return make_shared<AnalyzeTRAJ>(pair, args);
}

}  // namespace feasst

