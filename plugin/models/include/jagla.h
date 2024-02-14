
#ifndef FEASST_MODELS_JAGLA_H_
#define FEASST_MODELS_JAGLA_H_

#include <memory>
#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The Jagla model as described in
  https://doi.org/10.1103/PhysRevE.98.012138
  but with the variables renamed as follows:
  \f$\epsilon_1 \rightarrow \gamma\f$,
  \f$\epsilon_2 \rightarrow \epsilon\f$,
  \f$\lambda_1 \rightarrow \lambda\f$,
  \f$\lambda_2 \rightarrow r_c\f$.

  \f$ U(r) = \left\{
    \begin{array}{ll}
      \infty & : r < \sigma, \\
      \frac{\gamma(\lambda - r)-\epsilon(r-\sigma)}{\lambda-\sigma} & : \sigma \le r \le \lambda, \\
      -\frac{\epsilon(r_c - r)}{r_c-\lambda} & : \lambda \le r \le r_c, \\
      0 & : r \ge r_c.
    \end{array}
  \right. \f$
 */
class Jagla : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - num_discretized_steps: convert the continuous ramp potential into a number
      of dicontinuous, discretized steps.
      See Fig. 1(d) of https://doi.org/10.1103/PhysRevE.98.012138 .
      If 0, then use the linear ramp potential (default: 0).
   */
  explicit Jagla(argtype args = argtype());
  explicit Jagla(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the argument num_discretized_steps in read-only fashion.
  double num_discretized_steps() const { return num_discretized_steps_; }

  void precompute(const ModelParams& existing) override;

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<Jagla>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<Jagla>(args); }
  explicit Jagla(std::istream& istr);
  virtual ~Jagla() {}

  //@}
 private:
  double num_discretized_steps_;
  int lambda_index_ = -1;
  int gamma_index_ = -1;
};

inline std::shared_ptr<Jagla> MakeJagla(argtype args = argtype()) {
  return std::make_shared<Jagla>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_JAGLA_H_
