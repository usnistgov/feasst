
#ifndef FEASST_MONTE_CARLO_PERTURB_VOLUME_H_
#define FEASST_MONTE_CARLO_PERTURB_VOLUME_H_

#include "monte_carlo/include/perturb.h"

namespace feasst {

/**
  Change the volume of the system.
 */
class PerturbVolume : public Perturb {
 public:
  explicit PerturbVolume(argtype args = argtype());
  explicit PerturbVolume(argtype * args);

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false) override;

  /// Change volume
  void change_volume(const double delta_volume,
      System * system,
      const Select& select) {
    system->change_volume(delta_volume, args_); }

  void revert(System * system) override;
  void finalize(System * system) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbVolume(std::istream& istr);
  virtual ~PerturbVolume() {}

 private:
  argtype args_;

  // temporary
  double volume_change_;
};

inline std::shared_ptr<PerturbVolume> MakePerturbVolume(
    argtype args = argtype()) {
  return std::make_shared<PerturbVolume>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_VOLUME_H_
