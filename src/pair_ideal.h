/**
 * \file
 *
 * \brief ideal gas, no interactions
 *
 */
#ifndef PAIR_IDEAL_H_
#define PAIR_IDEAL_H_

#include "./pair.h"

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

  double peTot() { return 0.; }   //!< total potential energy of system
  double vrTot() { return 0.; }   //!< total virial of system

 protected:
};

#endif  // PAIR_IDEAL_H_

