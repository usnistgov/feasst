/**
 * Defines planar barriers which lie in a plane along two coordinate axes.
 */

#ifndef BARRIER_PLANAR_H_
#define BARRIER_PLANAR_H_

#include "./barrier.h"

namespace feasst {

class BarrierPlanar : public Barrier {
 public:
  /**
   * orthogonal planar barriers are defined by one coordinate point 
   * and direction, e.g., coord_ = 5, dimension_ = 0, and direction_ = 1
   * for wall at x >= 5
   */
  BarrierPlanar(const double coordinate, const int direction,
    const int dimen);
  virtual ~BarrierPlanar() {}

  /// \return the potential for a spherical object located at given coordinates
  double potential(const vector<double> coordinate,
    const double diameter);
  
 protected:
  double coord_;
  int direction_;
  int dimen_;
};

}  // namespace feasst

#endif  // BARRIER_PLANAR_H_

