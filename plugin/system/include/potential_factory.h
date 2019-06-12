
#ifndef FEASST_SYSTEM_POTENTIAL_FACTORY_H_
#define FEASST_SYSTEM_POTENTIAL_FACTORY_H_

#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include "system/include/potential.h"
//#include "utils/include/timer.h"

namespace feasst {

/**
 */
class PotentialFactory {
 public:
  PotentialFactory() {}

  void add_potential(const Potential potential) {
    potentials_.push_back(potential);
  }

  double energy(Configuration * config) {
    double en = 0;
    int index = 0;
    while ((index < static_cast<int>(potentials_.size())) and
           (en < NEAR_INFINITY)) {
      en += potentials_[index].energy(config);
      ++index;
    }
    DEBUG("en " << en);
    DEBUG(str());
    return en;
  }

  double energy(const Select& select, Configuration * config) {
    double en = 0;
    int index = 0;
    while ((index < static_cast<int>(potentials_.size())) and
           (en < NEAR_INFINITY)) {
      en += potentials_[index].energy(select, config);
      ++index;
    }
    DEBUG("en " << en);
    DEBUG(str());
    return en;
  }

  void revert() {
    for (Potential& potential : potentials_) {
      potential.revert();
    }
  }

  const std::vector<Potential>& potentials() const { return potentials_; }

  std::vector<double> stored_energy_profile() const {
    std::vector<double> en;
    for (const Potential& potential : potentials_) {
      en.push_back(potential.stored_energy());
    }
    return en;
  }

  double stored_energy() const {
    std::vector<double> en = stored_energy_profile();
    return std::accumulate(en.begin(), en.end(), 0.);
  }

  void precompute(Configuration * config) {
    for (Potential& potential : potentials_) {
      potential.precompute(config);
    }
  }

  int num() const { return static_cast<int>(potentials_.size()); }

  std::string str() const {
    std::stringstream ss;
    ss << "PotentialFactory: ";
    for (const Potential& potential : potentials_) {
      ss << potential.stored_energy() << " ";
    }
    return ss.str();
  }

  void serialize(std::ostream& sstr) const {
    feasst_serialize_version(1, sstr);
    feasst_serialize_fstobj(potentials_, sstr);
  }

  PotentialFactory(std::istream& sstr) {
    feasst_deserialize_version(sstr);
    feasst_deserialize_fstobj(&potentials_, sstr);
  }

 private:
  std::vector<Potential> potentials_;
//  Timer timer_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_POTENTIAL_FACTORY_H_
