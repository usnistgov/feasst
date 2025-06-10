
#ifndef FEASST_MODELS_F3C_H_
#define FEASST_MODELS_F3C_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
\rst
This class attempts to implement the F3C water model potential described in
:footcite:t:`levitt_calibration_1997`.
Unfortunately, we have not been able to reproduce those results, and there may
be an error in our implementation.
Use this implementation at your own risk.

The charge conversion factor assuming the following units:
  1. length: Angstroms
  2. energy: kJ/mol
  3. charge: elementary

When distance is near zero, returns a large positive number.

Acknowledgements to Samiha Sharlin for help with testing F3C.

References:

.. footbibliography::
\endrst
 */
class F3C : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - Asc: scaling factor (default: 1).
   */
  explicit F3C(argtype args = argtype());
  explicit F3C(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void precompute(const Configuration& config) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<F3C>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<F3C>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit F3C(std::istream& istr);
  virtual ~F3C() {}

  //@}
 private:
  double Asc_;
  double conversion_factor_ = 0.;
};

}  // namespace feasst

#endif  // FEASST_MODELS_F3C_H_
