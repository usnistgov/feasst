
#ifndef FEASST_CHAIN_RECENTER_PARTICLES_H_
#define FEASST_CHAIN_RECENTER_PARTICLES_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Periodically set particle positions to their geometric centers.
 */
class RecenterParticles : public ModifyUpdateOnly {
 public:
  /**
    args:
    - group_index: index of group to recenter (default: 0).
   */
  explicit RecenterParticles(const argtype &args = argtype());

  int group_index() const { return group_index_; }

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override {
    system->get_configuration()->recenter_particle_positions(group_index());
  }

  // serialize
  std::string class_name() const override {
    return std::string("RecenterParticles"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<RecenterParticles>(istr); }
  RecenterParticles(std::istream& istr);

 private:
  int group_index_;
};

inline std::shared_ptr<RecenterParticles> MakeRecenterParticles(
    const argtype &args = argtype()) {
  return std::make_shared<RecenterParticles>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_RECENTER_PARTICLES_H_
