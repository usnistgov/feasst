
#ifndef FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
#define FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Translate the positions of the selection uniformly random within postivie
  and negative tunable parameter values.
  The tunable parameter must be within the bounds of 2*NEAR_ZERO and half
  of Domain::max_side_length.
 */
class PerturbTranslate : public PerturbMove {
 public:
  //@{
  /** @name Arguments
    - dimension: if == -1, uniform translation in all dimensions (default: -1).
      Otherwise, translate only specifically +/- value (without tuning) in
      specific dimension.
   */
  explicit PerturbTranslate(argtype args = argtype());
  explicit PerturbTranslate(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Initialize minimum and maximum tunable parameter based on domain.
  void precompute(TrialSelect * select, System * system) override;

  /// Set the anchor to the current position of selection.
  void mid_stage(const TrialSelect& select, const System& system) override;
  void begin_stage(const TrialSelect& select) override;

  /// Change the position in the selection given a trajectory.
  void update_selection(const Position& trajectory,
      TrialSelect * select);

  /// Move the selected particles given a trajectory.
  void move(
      const Position& trajectory,
      System * system,
      TrialSelect * select);

  /// Move the selected particles.
  /// The particles are translated by +/- a maximum of the Tunable parameter.
  void move(const bool is_position_held, System * system, TrialSelect * select,
            Random * random, Acceptance * acceptance) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbTranslate(std::istream& istr);
  virtual ~PerturbTranslate() {}

  //@}
 private:
  int dimension_;

  // temporary objects not serialized.
  Position trajectory_, anchor_, new_pos_;
  bool anchor_set_ = false;
};

inline std::shared_ptr<PerturbTranslate> MakePerturbTranslate(
    argtype args = argtype()) {
  return std::make_shared<PerturbTranslate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
