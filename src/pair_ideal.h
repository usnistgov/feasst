#ifndef PAIR_IDEAL_H_
#define PAIR_IDEAL_H_

#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Ideal gas with no interactions.
 * PairIdeal still keeps track of neighbors for debugging purposes.
 */
class PairIdeal : public Pair {
 public:
  /// Constructor
  /// @param rCut HWH:depreciated. This variable has no use in this class.
  PairIdeal(Space* space, const double rCut);

  PairIdeal(Space* space, const char* fileName) : Pair(space, 0.) {
    ASSERT(0, "no restart implemented"); }
  ~PairIdeal();
  virtual PairIdeal* clone(Space* space) const {
    PairIdeal* p = new PairIdeal(*this); p->reconstruct(space); return p;
  }

  void initEnergy();     //!< function to calculate forces, given positions

  /// Potential energy of multiple particles.
  double multiPartEner(const vector<int> multiPart, const int flag);

 protected:
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_IDEAL_H_

