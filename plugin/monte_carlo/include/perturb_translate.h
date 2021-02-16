
#ifndef FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
#define FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Translate the positions of the selection.
 */
class PerturbTranslate : public PerturbMove {
 public:
  explicit PerturbTranslate(argtype args = argtype());
  explicit PerturbTranslate(argtype * args);

  /// Initialize minimum and maximum tunable parameter based on domain.
  void precompute(TrialSelect * select, System * system) override;

  /// Change the position in the selection given a trajectory.
  void update_selection(const Position& trajectory,
      TrialSelect * select);

  /// Move the selected particles given a trajectory.
  void move(
      const Position& trajectory,
      System * system,
      TrialSelect * select);

  /// Move the selected particles using the tuning parameter.
  void move(System * system, TrialSelect * select, Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbTranslate(std::istream& istr);
  virtual ~PerturbTranslate() {}

 private:
  // temporary objects
  Position trajectory_;
};

inline std::shared_ptr<PerturbTranslate> MakePerturbTranslate(
    argtype args = argtype()) {
  return std::make_shared<PerturbTranslate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
