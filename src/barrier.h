/**
 * This class serves as a container for derives classes with specific barrier
 * implementations, in order to interface with PairWall.
 */

#ifndef BARRIER_H_
#define BARRIER_H_

#include "./base.h"

namespace feasst {

class Barrier : public Base {
 public:
  Barrier();    //!< Constructor
  virtual ~Barrier() {}

  /// \return the potential for a spherical object located at given coordinates
  virtual double potential(const vector<double> coordinate,
    const double diameter);
  
  /// add a barrier
  void add(shared_ptr<Barrier> barrier) { 
    barriers_.push_back(barrier);
  }

  /**
   * Add a barrier which lies in a plane along two coordinate axes.
   * The other axes, given by dimen, has a coordinate and direction
   * where the barrier is present.
   */
  void addOrthogonalPlanar(const double coord, const int direction,
    const int dimen);

 protected:
  vector<shared_ptr<Barrier> > barriers_;
};

}  // namespace feasst

#endif  // BARRIER_H_

