
#ifndef FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Rigid translation of a cluster of particles.
 */
class TrialTranslateCluster : public Trial {
 public:
  /// These arguments are sent to both PerturbTranslate and TrialStage.
  TrialTranslateCluster(std::shared_ptr<NeighborCriteria> neighbor_criteria,
                        const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialTranslateCluster(std::istream& istr);
  virtual ~TrialTranslateCluster() {}

 protected:
  void serialize_trial_translate_cluster_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialTranslateCluster> MakeTrialTranslateCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialTranslateCluster>(neighbor_criteria, args);
}

/**
  Rigid translation of a cluster of particles.
 */
class TrialRotateCluster : public Trial {
 public:
  /// These arguments are sent to both PerturbRotate and TrialStage.
  TrialRotateCluster(std::shared_ptr<NeighborCriteria> neighbor_criteria,
                     const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialRotateCluster(std::istream& istr);
  virtual ~TrialRotateCluster() {}

 protected:
  void serialize_trial_rotate_cluster_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialRotateCluster> MakeTrialRotateCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialRotateCluster>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_
