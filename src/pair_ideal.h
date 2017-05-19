/**
 * \file
 *
 * \brief ideal gas, no interactions
 *
 */
#ifndef PAIR_IDEAL_H_
#define PAIR_IDEAL_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairIdeal : public Pair {
 public:
  PairIdeal(Space* space, const double rCut);
  ~PairIdeal();
  virtual PairIdeal* clone(Space* space) const {
    PairIdeal* p = new PairIdeal(*this); p->reconstruct(space); return p;
  }

  int initEnergy();     //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

 protected:
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_IDEAL_H_

