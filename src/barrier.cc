#include "./barrier.h"
#include "./barrier_planar.h"

namespace feasst {

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

}  // namespace feasst

