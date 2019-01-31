
#ifndef FEASST_CORE_MODEL_H_
#define FEASST_CORE_MODEL_H_

#include "core/include/particle.h"
#include "core/include/configuration.h"

namespace feasst {

class VisitModel;

class Model {
 public:
  /// Visit the model over the entire configuration.
  /// Optionally, restrict to groups of given index.
  virtual double compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) const = 0;

  /// Same as above, except the group index is assumed to be zero (which is all
  /// particles and sites in the configuration).
  virtual double compute(
    Configuration * config,
    VisitModel * visitor) const = 0;

  /// Visit the model over a selection of the configuration.
  /// Optionally, restrict to groups of given index, which is only relevant for
  /// multibody models (e.g., two body and not one body).
  virtual double compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) const = 0;

  /// Same as above, except the group index is assumed to be zero (which is all
  /// particles and sites in the configuration)
  virtual double compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) const = 0;
  virtual ~Model() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_H_
