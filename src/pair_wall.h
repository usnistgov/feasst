/**
 * Compute interactions of particles with walls defined by the Barrier class.
 */

#ifndef PAIR_WALL_H_
#define PAIR_WALL_H_

#include "./pair.h"
#include "./barrier.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairWall : public Pair {
 public:
  PairWall(Space* space, Barrier* barrier);
  ~PairWall() {}
  virtual PairWall* clone(Space* space) const {
    PairWall* p = new PairWall(*this); p->reconstruct(space); return p;
  }

  /*
   * \return potential energy of a selection of particles multiPart, and store
   * this the Pair class variable peSRone
   */
  double multiPartEner(const vector<int> multiPart, const int flag = 0);

  /*
   * Compute the potential energy of all particles and store this in the Pair
   * class variable peTot
   */
  int initEnergy();

  /// stores, restores or updates variables to avoid order recompute
  //   of entire configuration after every change
  virtual void update(const vector<int> mpart, const int flag,
                      const char* uptype);

 protected:
  Barrier* barrier_;
  double peWall_;
  double deWall_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_WALL_H_

