
#ifndef FEASST_CORE_MODEL_H_
#define FEASST_CORE_MODEL_H_

#include "core/include/particle.h"
#include "core/include/configuration.h"

namespace feasst {

class VisitModel;

class Model {
 public:
  /// Visit the model over the entire configuration.
  virtual double compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) const = 0;
  virtual double compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) const = 0;

  /// Same as above, except the group index is assumed to be zero (which is all
  /// particles and sites in the configuration).
  virtual double compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) const = 0;
  virtual double compute(
    Configuration * config,
    VisitModel * visitor) const = 0;

  /// Visit the model over a selection of the configuration.
  /// Optionally, restrict to groups of given index, which is only relevant for
  /// multibody models (e.g., two body and not one body).
  virtual double compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) const = 0;
  virtual double compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) const = 0;

  /// Same as above, except the group index is assumed to be zero (which is all
  /// particles and sites in the configuration)
  virtual double compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) const = 0;
  virtual double compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) const = 0;
  virtual ~Model() {}

  /// Precompute model parameters based on existing model parameters.
  virtual void precompute(const ModelParams& existing) {}

  // https://isocpp.org/wiki/faq/serialization
  //typedef std::shared_ptr<Model> (*Factory)(std::istream&);

  /// Output a serialized version of the existing model.
  virtual void serialize(std::ostream& ostr) const { ERROR("not implemented"); }

  // Derived class implementation of a serialization.
  virtual std::shared_ptr<Model> create(std::istream& istr) const {
    ERROR("not implemented"); }

  // Returns a static mapping of class name to model.
  std::map<std::string, std::shared_ptr<Model> >& deserialize_map();

  /// Return a model given a serialization.
  std::shared_ptr<Model> deserialize(std::istream& istr) {
    return template_deserialize(deserialize_map(), istr);
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_H_
