
#ifndef FEASST_MONTE_CARLO_PERTURB_ANYWHERE_H_
#define FEASST_MONTE_CARLO_PERTURB_ANYWHERE_H_

#include <memory>
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/// Rigidly move anywhere in the box with any orientation.
class PerturbAnywhere : public PerturbMove {
 public:
  PerturbAnywhere();

  /// Set the particle position of select to center.
  void set_position(const Position& center,
                    System * system,
                    TrialSelect * select);

  void move(System * system, TrialSelect * select, Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbAnywhere(std::istream& istr);
  virtual ~PerturbAnywhere() {}

 private:
  PerturbRotate rotate_;

  // temporary
  PerturbTranslate translate_;
  Position random_in_box_;
};

inline std::shared_ptr<PerturbAnywhere> MakePerturbAnywhere() {
  return std::make_shared<PerturbAnywhere>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ANYWHERE_H_
