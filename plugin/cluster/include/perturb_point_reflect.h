
#ifndef FEASST_CLUSTER_PERTURB_POINT_REFLECT_H_
#define FEASST_CLUSTER_PERTURB_POINT_REFLECT_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Translate the positions of the selection.
 */
class PerturbPointReflect : public PerturbMove {
 public:
  explicit PerturbPointReflect(const argtype& args = argtype());

  /// Initialize minimum and maximum tunable parameter based on domain.
  void precompute(TrialSelect * select, System * system) override;

  /// Change the position in the selection given a point reflection.
  void update_selection(const Position& reflect,
      TrialSelect * select);

  /// Move the selected particles given a point reflection.
  void move(
      const Position& reflect,
      System * system,
      TrialSelect * select);

  /// Move the selected particles using the tuning parameter.
  void move(System * system, TrialSelect * select, Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbPointReflect(std::istream& istr);
  virtual ~PerturbPointReflect() {}

 private:
  // temporary objects
  Position reflect_;
};

inline std::shared_ptr<PerturbPointReflect> MakePerturbPointReflect(
    const argtype& args = argtype()) {
  return std::make_shared<PerturbPointReflect>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_PERTURB_POINT_REFLECT_H_
