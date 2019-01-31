
#ifndef FEASST_CORE_PERTURB_CONFIGS_H_
#define FEASST_CORE_PERTURB_CONFIGS_H_

#include "core/include/perturb.h"

namespace feasst {

/**
  Originally implemented for testing purpose and limited scope.
 */
class PerturbConfigs : public Perturb {
 public:
  void before_attempt() override {
    Perturb::before_attempt();
  }

  const Select& selection() const override {
    ERROR("not implemented");
    return donor_part_;
  }

  void transfer_particle(const int index, System * system, const int donor, const int acceptor) {
    ASSERT(donor != acceptor, "donor and acceptor cannot be the same");
    donor_part_.particle(index, system->configuration(donor));
    Configuration * donor_config_ = system->get_configuration(donor);
    Configuration * acceptor_config_ = system->get_configuration(acceptor);
    const Particle& part = donor_part_.particle(system->configuration(donor));
    const int type = part.type();
    acceptor_config_->add_particle(type);
    acceptor_part_.last_particle_added(acceptor_config_);
    acceptor_config_->replace_position(acceptor_part_, part);
    donor_config_->remove_particle(donor_part_);
  }

  void revert() override { ERROR("not implemented"); }

  ~PerturbConfigs() {}

 private:
  SelectList donor_part_;
  SelectList acceptor_part_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_CONFIGS_H_
