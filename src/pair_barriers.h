/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 * Compute interactions of particles with walls defined by the Barrier class.
 */

#ifndef PAIR_BARRIERS_H_
#define PAIR_BARRIERS_H_

#include "./pair.h"
#include "./barrier.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairBarriers : public Pair {
 public:
  PairBarriers (Space* space, SpeciesBarriers* specbarriers) : Pair(space, 0.) {}
  ~PairBarriers () {}
  virtual PairBarriers* clone (Space* space) const {
    PairBarriers* p = new PairBarriers (*this); p->reconstruct(space); return p;
  }

  /*
   * \return potential energy of a selection of particles multiPart, and store
   * this the Pair class variable peSRone
   */
  double multiPartEner(const vector<int> multiPart, const int flag) {
    if (flag == 0 || multiPart.size() < 2) {} return 0;
  }

  /*
   * Compute the potential energy of all particles and store this in the Pair
   * class variable peTot
   */
  virtual int initEnergy() { return 0; }
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_BARRIERS_H_
