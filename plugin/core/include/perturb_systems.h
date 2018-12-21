
#ifndef FEASST_CORE_PERTURB_SYSTEMS_H_
#define FEASST_CORE_PERTURB_SYSTEMS_H_

#include "core/include/perturb.h"

namespace feasst {

/**
  Originally implemented for testing purpose and limited scope.
 */
class PerturbSystems : public Perturb {
 public:
  void before_attempt() override {
    Perturb::before_attempt();
  }

  const Select& selection() const override {
    ERROR("not implemented");
    return donor_part_;
  }

  void transfer_particle(const int index, System * donor, System * acceptor) {
    donor_part_.particle(index, donor->config());
    Configuration * donor_config_ = donor->configuration(0);
    Configuration * acceptor_config_ = acceptor->configuration(0);
    const Particle& part = donor_part_.particle(*donor->configuration(0));
    const int type = part.type();
    acceptor_config_->add_particle(type);
    acceptor_part_.last_particle_added(acceptor_config_);
    acceptor_config_->replace_position(acceptor_part_, part);
    donor_config_->remove_particle(donor_part_);
  }

  void revert() { ERROR("not implemented"); }

  ~PerturbSystems() {}

 private:
  SelectList donor_part_;
  SelectList acceptor_part_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_SYSTEMS_H_
