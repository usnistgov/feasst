#include "./barrier.h"
#include "./barrier_planar.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Barrier::Barrier() {
  className_.assign("Barrier");
}

double Barrier::potential(const vector<double> coordinate,
    const double diameter) {
  double pe = 0.;
  for (int ibarr = 0; ibarr < static_cast<int>(barriers_.size()); ++ibarr) {
    pe += barriers_[ibarr]->potential(coordinate, diameter);
  }
  return pe;
}

void Barrier::addOrthogonalPlanar(const double coord, const int direction,
  const int dimen) {
  shared_ptr<BarrierPlanar> newBarrier = make_shared<BarrierPlanar>(
    coord, direction, dimen);
  add(newBarrier);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

