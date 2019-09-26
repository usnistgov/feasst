
#ifndef FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
#define FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Translate the positions of the selection.
 */
class PerturbTranslate : public PerturbMove {
 public:
  PerturbTranslate(const argtype& args = argtype()) : PerturbMove(args) {
    class_name_ = "PerturbTranslate";
  }

  void precompute(TrialSelect * select, System * system) override {
    set_tunable_min_and_max(2*NEAR_ZERO,
      0.5*system->configuration().domain().max_side_length());
  }

  /// Change the position in the selection given a trajectory.
  void update_selection(const Position& trajectory,
      TrialSelect * select) {
    SelectList * displaced = select->get_mobile();
    for (int select_index = 0;
         select_index < displaced->num_particles();
         ++select_index) {
      Position displaced_part(displaced->particle_positions()[select_index]);
      displaced_part.add(trajectory);
      displaced->set_particle_position(select_index, displaced_part);
      for (int site = 0;
           site < static_cast<int>(displaced->site_indices(select_index).size());
           ++site) {
        Position displaced_site(displaced->site_positions()[select_index][site]);
        displaced_site.add(trajectory);
        displaced->set_site_position(select_index, site, displaced_site);
      }
    }
  }

  /// Move the selected particles given a trajectory.
  void move(
      const Position& trajectory,
      System * system,
      TrialSelect * select) {
    update_selection(trajectory, select);
    system->get_configuration()->update_positions(select->mobile());
  }

  /// Move the selected particles using the tuning parameter.
  void move(System * system, TrialSelect * select, Random * random) override {
    random->position_in_cube(
      system->dimension(),
      tunable().value(),
      &trajectory_
    );
    DEBUG("max move " << tunable().value());
    ASSERT(tunable().value() > NEAR_ZERO, "tunable(" << tunable().value()
      << ") is too small");
    move(trajectory_, system, select);
  }

  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbTranslate(std::istream& istr);
  virtual ~PerturbTranslate() {}

 private:
  // temporary objects
  Position trajectory_;
};

inline std::shared_ptr<PerturbTranslate> MakePerturbTranslate(const argtype& args = argtype()) {
  return std::make_shared<PerturbTranslate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
