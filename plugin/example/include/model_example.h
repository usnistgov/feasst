
#ifndef FEASST_EXAMPLE_MODEL_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_EXAMPLE_H_

#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  Users may implement their own custom Model using this template file.
  This file is located at /path/to/feasst/plugin/example/include/model_example.h
  In C++, the header file has the extension of ".h" and is typically used to
  describe the interface.
  The implementation is found in
  /path/to/feasst/plugin/example/src/model_example.cpp
  Both h and cpp files correspond to one C++ object.

  To implement your custom ModelTwoBody (pair-wise Model), in most cases you
  will simply need to change the implementation of the energy(...) function in
  the cpp file.

  The inputs to ModelExample::energy(...) are more difficult to change, because
  they are set by ModelTwoBody, which is the parent class of ModelExample.
  ModelExample simply overrides the energy(...) function.

  Arguments can however be added to the ModelExample constructor, as shown
  below.

  To rename ModelExample, which is also necessary for implementing multiple
  different models, you need to do two things.
  1. Replace all text ModelExample with the new name in both the cpp and h file.
     With regex in VIM, this can be done as :%s/ModelExample/NewName/g
  2. Replace all text FEASST_EXAMPLE_MODEL_EXAMPLE_H_ with
     FEASST_EXAMPLE_NEW_NAME_H_ in order to update the header guards.
  With these two steps, you can now have multiple objects or pairs of h and cpp
  files.
 */
class ModelExample : public ModelTwoBody {
 public:
  /**
    Use this comment to describe any user arguments (or delete if no user
    arguments are required).

    args:
    - example_argument: [put description of an example argument]
      (default: 0).

    Remember to put default value here that is used if no argument is provided.
    Make sure the default value corresponds to the one in the .cpp file.
   */
  explicit ModelExample(argtype args = argtype());
  explicit ModelExample(argtype * args);

  /// Return the example_argument in a constant, read-only fashion.
  double example_argument() const { return example_argument_; }

  /**
    The energy between two site types depends upon the distance between
    them and the model parameters.
    Find other existing ModelTwoBody, such as LennardJones or
    LennardJonesCutShift for more inspiration.
    Then, implement your own model in /path/to/feasst/plugin/example/src/model_example.cpp
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
  double example_argument_;
};

inline std::shared_ptr<ModelExample> MakeModelExample(argtype args = argtype()) {
  return std::make_shared<ModelExample>(args);
}

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_EXAMPLE_H_
