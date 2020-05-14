/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef ANALYZE_TRAJ_H_
#define ANALYZE_TRAJ_H_

#include "./analyze.h"

namespace feasst {

/**
 * Compute cluster properties
 */
class AnalyzeTRAJ : public Analyze {
 public:
  /// Constructor
  AnalyzeTRAJ(Pair *pair,
    /**
     * allowed string key pairs (e.g., dictionary):
     *
     *  format : trajectory format (default: xyz).
     *   - xyz : see Space::readXYZ and Space::printXYZ
     */
    const argtype &args = argtype());

  /// Write file.
  void write();

  // Write file the same way, regardless of macrostate
  void write(CriteriaWLTMMC *c) { write(); }

  /// Initialize production flag. 1 is on, 0 is off. Default is 0.
  void initProduction(const int flag = 0) { if (flag == 1) firstFlag_ = 1; }

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  AnalyzeTRAJ(Pair *pair, const char* fileName);

  ~AnalyzeTRAJ() {}
  AnalyzeTRAJ* clone(Pair* pair) const {
    AnalyzeTRAJ* a = new AnalyzeTRAJ(*this);
    a->reconstruct(pair); return a;
  }
  shared_ptr<AnalyzeTRAJ> cloneShrPtr(Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeTRAJ, Analyze>(cloneImpl(pair)));
  }

 protected:
  void defaultConstruction_();
  int firstFlag_ = 1;   // set to zero after first write
  std::string format_;

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Pair *pair) const {
    shared_ptr<AnalyzeTRAJ> a = make_shared<AnalyzeTRAJ>(*this);
    a->reconstruct(pair); return a;
  }
};

/// Factory method
shared_ptr<AnalyzeTRAJ> makeAnalyzeTRAJ(Pair *pair,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // ANALYZE_TRAJ_H_

