
#ifndef FEASST_SYSTEM_MODEL_H_
#define FEASST_SYSTEM_MODEL_H_

#include <memory>
#include <string>
#include <map>

namespace feasst {

typedef std::map<std::string, std::string> argtype;
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

  /// Return the number of bodies in the model (e.g., two-body/pairwise = 2).
  virtual int num_body() const = 0;

  /// Precompute model parameters based on existing model parameters.
  virtual void precompute(const ModelParams& existing);

  /// Return the ModelParams index of epsilon.
  int epsilon_index() const { return epsilon_index_; }

  /// Return the ModelParams index of sigma.
  int sigma_index() const { return sigma_index_; }

  /// Return the ModelParams index of cutoff.
  int cutoff_index() const { return cutoff_index_; }

  /// Return the ModelParams index of charge.
  int charge_index() const { return charge_index_; }

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

  // Derived class implementation of a serialization.
  virtual std::shared_ptr<Model> create(argtype * args) const;

  // Returns a static mapping of class name to model.
  std::map<std::string, std::shared_ptr<Model> >& deserialize_map();

  /// Return a model given a serialization.
  std::shared_ptr<Model> deserialize(std::istream& istr);

  /// Return a model given arguments.
  std::shared_ptr<Model> factory(const std::string name, argtype * args);

  // Constructor reads class_name
  explicit Model(std::istream& istr);

 protected:
  std::string class_name_;
  void serialize_model_(std::ostream& ostr) const;

 private:
  int epsilon_index_ = -1;
  int sigma_index_ = -1;
  int cutoff_index_ = -1;
  int charge_index_ = -1;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_H_
