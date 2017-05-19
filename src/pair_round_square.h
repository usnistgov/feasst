/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef PAIR_ROUND_SQUARE_H_
#define PAIR_ROUND_SQUARE_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairRoundSquare : public Pair {
 public:
  PairRoundSquare(Space* space, const double rCut);
  PairRoundSquare(Space* space, const char* fileName);
  ~PairRoundSquare() {}
  virtual PairRoundSquare* clone(Space* space) const {
    PairRoundSquare* p = new PairRoundSquare(*this); p->reconstruct(space); return p;
  }

  int initEnergy() { return 0; }

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag) {
    return multiPart.size()*flag;
  }

 protected:
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_ROUND_SQUARE_H_

