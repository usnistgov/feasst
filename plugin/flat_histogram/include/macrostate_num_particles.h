
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the total number of particles in the system.
 */
class MacrostateNumParticles : public Macrostate {
 public:
  MacrostateNumParticles(const Histogram& histogram,
    const argtype& args = argtype()) : Macrostate(histogram, args) {}

  double value(const System* system, const Criteria* criteria) override {
    return static_cast<double>(system->configuration().num_particles());
  }

  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateNumParticles(std::istream& istr);
  virtual ~MacrostateNumParticles() {}

 private:
  const std::string class_name_ = "MacrostateNumParticles";
};

inline std::shared_ptr<MacrostateNumParticles> MakeMacrostateNumParticles(
    const Histogram& histogram, const argtype& args = argtype()) {
  return std::make_shared<MacrostateNumParticles>(histogram, args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_
