
#ifndef FEASST_EXAMPLE_MODEL_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_EXAMPLE_H_

#include <memory>
#include <map>
#include <string>
#include "system/include/model_two_body.h"

namespace feasst {

class ModelParams;

typedef std::map<std::string, std::string> argtype;

// Comments such as this beginning with // or /* do not appear in the
// HTML/PDF documentation.
// Note that references such as "jagla_core-softened_1999" require the rst
// environment and an entry in /feasst/dev/sphinx/refs.bib
// However, rst environment (beginning \rst ending \endrst) does not
// link class names.
// For this reason, we demonstrate here a mixture of documentation in both the
// default doxygen and, further down, the \rst environment.
// Comments beginning with /// or /** appear in documentation as follows.

/**
  Add a model to FEASST by using this file as a template and instruction set.
  Follow the same steps detailed in /feasst/plugin/example/README.rst.
  In summary, copy model_example.[h/cpp] to new_name.[h/cpp], replace
  ModelExample with NewName, then replace MODEL_EXAMPLE with NEW_NAME.

  In some cases, TablePotential may be a more convenient and efficient approach
  to implementing a custom potential in FEASST.

  For more inspiration, take a look at other existing ModelTwoBody, such as
  LennardJones or LennardJonesCutShift, or the models plugin.
  Sometimes it is easiest to find an existing model that is the more similar to
  what you need, and copy/rename those files instead.

\rst
As an example, consider the Jagla model,\ :footcite:p:`jagla_core-softened_1999`
but with the variables renamed as follows:
:math:`\epsilon_1 \rightarrow \gamma`,
:math:`\epsilon_2 \rightarrow \epsilon`,
:math:`\lambda_1 \rightarrow \lambda`,
:math:`\lambda_2 \rightarrow r_c`.

.. math::

    U(r) = \left\{
      \begin{array}{ll}
        \infty & : r < \sigma, \\
        \frac{\gamma(\lambda - r)-\epsilon(r-\sigma)}{\lambda-\sigma} & : \sigma \le r \le \lambda, \\
        -\frac{\epsilon(r_c - r)}{r_c-\lambda} & : \lambda \le r \le r_c, \\
        0 & : r \ge r_c.
      \end{array}
    \right.

In some studies, the continuous linear ramp has been converted into a number of
discrete steps to enable simulations with discontinuous molecular dynamics.
For example, see Fig. 1(d) of :footcite:t:`de_haro_structural_2018`

References:

.. footbibliography::
\endrst*/
class ModelExample : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - num_discretized_steps: convert the continuous ramp potential into a number
      of dicontinuous, discretized steps.
      If 0, then use the linear ramp potential (default: 0).

    Use this comment to describe and add any user arguments (or delete if no
    user arguments are required).
    Ensure the default value corresponds to the one in the .cpp file.
    Note that "Site Properties" in fstprt files are not added here.
   */
  explicit ModelExample(argtype args = argtype());
  explicit ModelExample(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the argument num_discretized_steps in read-only fashion.
  double num_discretized_steps() const { return num_discretized_steps_; }

  /*
    Precompute is called before the simulation starts and is used to update
    the ModelParams, including their mixing, which could include relatively
    costly sqrt operations.
    Precompute is also searches the names of the model parameters for their
    appropriate numerical index, so they can be more quickly accessed by
    index during the simulation.
   */
  void precompute(const Configuration& config) override;

  /**
    The energy between two site types depends upon the distance between
    them and the model parameters.
   */
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  /// Serialization is the flatening of an object into a stream of characters
  /// which may be saved to file and later deserialized back into the object.
  void serialize(std::ostream& ostr) const override;

  /// Deserialization is the reverse process of the above.
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelExample>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelExample>(args); }
  explicit ModelExample(std::istream& istr);
  virtual ~ModelExample() {}

  //@}
 private:
  /**
    Declare private member variables here.
    Private member variables are only accessible within the scope of this class.
   */
  double num_discretized_steps_;
  int lambda_index_ = -1;
  int gamma_index_ = -1;
};

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_EXAMPLE_H_
