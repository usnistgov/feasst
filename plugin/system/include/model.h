
#ifndef FEASST_SYSTEM_MODEL_H_
#define FEASST_SYSTEM_MODEL_H_

#include <memory>
#include <string>
#include <map>

namespace feasst {

class ModelParams;
class Select;
class Configuration;
class VisitModel;

class Model {
 public:
  Model() { class_name_ = "Model"; }

  /// Visit the model over the entire configuration.
  virtual double compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) = 0;
  virtual double compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) = 0;

  /// Same as above, except the group index is assumed to be zero (which is all
  /// particles and sites in the configuration).
  virtual double compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) = 0;
  virtual double compute(
    Configuration * config,
    VisitModel * visitor) = 0;

  /// Visit the model over a selection of the configuration.
  /// Optionally, restrict to groups of given index, which is only relevant for
  /// multibody models (e.g., two body and not one body).
  virtual double compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) = 0;
  virtual double compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) = 0;

  /// Same as above, except the group index is assumed to be zero (which is all
  /// particles and sites in the configuration)
  virtual double compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) = 0;
  virtual double compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) = 0;
  virtual ~Model() {}

  /// Precompute model parameters based on existing model parameters.
  virtual void precompute(const ModelParams& existing) {}

  // Moved from ModelTwoBody to Model for ease of serialization of
  // ModelTwoBodyFactory
  virtual double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params);

  const std::string& class_name() const { return class_name_; }

  /// Output a serialized version of the existing model.
  virtual void serialize(std::ostream& ostr) const;

  // Derived class implementation of a serialization.
  virtual std::shared_ptr<Model> create(std::istream& istr) const;

  // Returns a static mapping of class name to model.
  std::map<std::string, std::shared_ptr<Model> >& deserialize_map();

  /// Return a model given a serialization.
  std::shared_ptr<Model> deserialize(std::istream& istr);

  // Constructor reads class_name
  explicit Model(std::istream& istr);

 protected:
  std::string class_name_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_H_
