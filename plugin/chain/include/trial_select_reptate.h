
#ifndef FEASST_CHAIN_TRIAL_SELECT_REPTATE_H_
#define FEASST_CHAIN_TRIAL_SELECT_REPTATE_H_

#include "chain/include/trial_select_end_segment.h"

namespace feasst {

const argtype TrialSelectReptateArg_ = {{"max_length", "1"}};

/// Select a random end point for reptation.
class TrialSelectReptate : public TrialSelectEndSegment {
 public:
  TrialSelectReptate() : TrialSelectEndSegment(TrialSelectReptateArg_) {
    class_name_ = "TrialSelectReptate";
  }

  void precompute(System * system) override {
    anchor_.clear();
    anchor_.add_site(0, 0);
    bonded_to_.clear();
    bonded_to_.add_site(0, 0);
  }

  void update_anchor(const bool is_endpoint_beginning,
    const System * system) override {
    const int particle_index = mobile_.particle_indices()[0];
    const Configuration& config = system->configuration();
    const Particle& particle = config.select_particle(particle_index);
    int anchor_index = 0;
    int site_bonded_to = particle.num_sites() - 2;
    DEBUG("is_endpoint_beginning " << is_endpoint_beginning);
    if (is_endpoint_beginning) {
      anchor_index = particle.num_sites() - 1;
      site_bonded_to = 1;
    }
    // for the old configuration, set the anchor to the old bond.
    anchor_.set_site(0, 0, anchor_index);
    anchor_.set_particle(0, particle_index);
    ASSERT(bonded_to_.replace_indices(particle_index, {site_bonded_to}),
      "bonded_to_ wasn't initialized to proper size on precompute");
  }

  void mid_stage() override {
    // exclude the anchor from interactions.
    // include interactions with site that use to be bonded
    mobile_.set_new_bond(anchor_);
    mobile_.set_old_bond(bonded_to_);
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectReptate(std::istream& istr);
  virtual ~TrialSelectReptate() {}

 protected:
  void serialize_trial_select_reptate_(std::ostream& ostr) const;

 private:
  // temporary
  SelectList bonded_to_;
};

inline std::shared_ptr<TrialSelectReptate> MakeTrialSelectReptate() {
  return std::make_shared<TrialSelectReptate>();
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_REPTATE_H_
