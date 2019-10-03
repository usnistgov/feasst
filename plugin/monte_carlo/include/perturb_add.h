
#ifndef FEASST_MONTE_CARLO_PERTURB_ADD_H_
#define FEASST_MONTE_CARLO_PERTURB_ADD_H_

#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

/**
  Add a particle to the system.
 */
// HWH optimize -> update cell list in finalize?
class PerturbAdd : public Perturb {
 public:
  PerturbAdd(const argtype& args = argtype());

  //initialize ghost selection in TrialSelect?
  void precompute(TrialSelect * select, System * system) override {
    select->set_ghost(true); }

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false
      ) override {
    add(system, select, random, empty_, is_position_held);
  }

  void add(
    System * system,
    TrialSelect * select,
    Random * random,
    /// place particle anywhere if center is of zero dimension.
    const Position& center,
    const bool is_position_held = false
  ) {
    DEBUG("is_position_held " << is_position_held);
    Configuration* config = system->get_configuration();
    config->revive(select->mobile());
    select->set_trial_state(2);

    // obtain probability
    const int particle_type = config->select_particle(
      select->mobile().particle_index(0)
    ).type();
    DEBUG("type " << particle_type);
    for (const Select& ghost : config->ghosts()) {
      DEBUG("ghost " << ghost.str());
    }
    select->set_probability(
      1./static_cast<double>(config->num_particles_of_type(particle_type)));

    if (center.dimension() == 0) {
      anywhere_.perturb(system, select, random, is_position_held);
    } else {
      anywhere_.set_position(center, system, select);
    }
    set_revert_possible(true, select);
  }

  void revert(System * system) override {
    DEBUG("revert_possible " << revert_possible());
    if (revert_possible()) {
      DEBUG(revert_select()->mobile().str());
      DEBUG("nump " << system->configuration().num_particles());
      system->get_configuration()->remove_particles(revert_select()->mobile());
      system->revert();
    }
  }

  std::string status_header() const override {
    std::stringstream ss;
    return ss.str();
  }

  std::string status() const override {
    std::stringstream ss;
    return ss.str();
  }

  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbAdd(std::istream& istr);
  virtual ~PerturbAdd() {}

 private:
  // temporary
  Position empty_;
  PerturbAnywhere anywhere_;
};

inline std::shared_ptr<PerturbAdd> MakePerturbAdd(const argtype& args = argtype()) {
  return std::make_shared<PerturbAdd>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ADD_H_
