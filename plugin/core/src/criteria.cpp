
#include "core/include/criteria.h"

namespace feasst {

void Criteria::set_beta(const double beta) {
  beta_ = beta;
  beta_initialized_ = true;
}

double Criteria::beta() const {
  ASSERT(beta_initialized_, "beta must be initialized before use");
  return beta_;
}

double Criteria::activity(const int particle_type) const {
  ASSERT(particle_type < static_cast<int>(activity_.size()),
    "activity of type(" << particle_type << ") must be initalized before use"); 
  return activity_[particle_type];
}

}  // namespace feasst
