
#ifndef FEASST_EXAMPLE_MODEL_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_EXAMPLE_H_

#include <memory>
#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"

namespace feasst {

/**
  Add a model to FEASST by using this file as a template and instruction set.

  1. Copy and paste this h file, and its corresponding cpp file with a new name:

     cd /path/to/feasst/plugin/example

     cp include/model_example.h include/new_name.h

     cp src/model_example.cpp src/new_name.cpp

  2. Replace ModelExample with NewName in the cpp and h files.
     With regex in VIM, this can be done as :%s/ModelExample/NewName/g

  3. Similarly, replace FEASST_EXAMPLE_MODEL_EXAMPLE_H_ with
     FEASST_EXAMPLE_NEW_NAME_H_ .

  4. Update the interface (feasst.h/feasst.i) and HTML documentation with the
     following BASH command:
     cd /path/to/feasst/build; python ../py/depend-py -s ../

  5. (optional) Add tests by performing steps 1 and 2 for
     /feasst/plugin/example/test/model_example.cpp .

  6. (optional) Change the serialization version found in two places in the cpp
     file from the current 5023 to something else (in both places).

  If an update to the code via git is required, note that step 4 will change the
  existing feasst.[i/h] files in the git repo.
  Use "git checkout" on feasst.[i/h] before "git pull" and then repeat step 4
  afterward.

  In some cases, TablePotential may be a more convenient and efficient approach
  to implementing a custom potential in FEASST.

  As an example, consider the Jagla model as described in
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

  For more inspiration, take a look at other existing ModelTwoBody, such as
  LennardJones or LennardJonesCutShift, or the models plugin.
  Sometimes it is easiest to find an existing model that is the more similar to
  what you need, and copy/rename those files instead.
 */
class ModelExample : public ModelTwoBody {
 public:
  /**
    args:
    - num_discretized_steps: convert the continuous ramp potential into a number
      of dicontinuous, discretized steps.
      See Fig. 1(d) of https://doi.org/10.1103/PhysRevE.98.012138 .
      If 0, then use the linear ramp potential (default: 0).

    Use this comment to describe and add any user arguments (or delete if no
    user arguments are required).
    Ensure the default value corresponds to the one in the .cpp file.
    Note that "Site Properties" in fstprt files are not added here.
   */
  explicit ModelExample(argtype args = argtype());
  explicit ModelExample(argtype * args);

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
  void precompute(const ModelParams& existing) override;

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

 private:
  /**
    Declare private member variables here.
    Private member variables are only accessible within the scope of this class.
   */
  double num_discretized_steps_;
  int lambda_index_ = -1;
  int gamma_index_ = -1;
};

inline std::shared_ptr<ModelExample> MakeModelExample(argtype args = argtype()) {
  return std::make_shared<ModelExample>(args);
}

/**
  While sigma, epsilon and cutoff are already implemented in Model,
  and the lambda parameter was already implemented in LennardJonesAlpha,
  the gamma parameter is not yet defined.
  Thus, we will go through the exercise here.
 */
class Gamma : public ModelParam {
 public:
  Gamma() : ModelParam() { class_name_ = "gamma"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Gamma>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit Gamma(std::istream& istr);
  virtual ~Gamma() {}
 private:
  double mix_(const double value1, const double value2) override;
};

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_EXAMPLE_H_
