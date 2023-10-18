
#ifndef FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt TrialTranslateCluster and TrialRotateCluster with equal probability.
 */
class TrialRigidCluster : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - rotate_param: initial value of the tunable parameter (default: 25).
    - translate_param: initial value of the tunable parameter (default: 0.1).
    - Trial arguments.
    - SelectCluster arguments.
    - Tunable arguments.
   */
  explicit TrialRigidCluster(argtype args = argtype());
  explicit TrialRigidCluster(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialRigidCluster>(args); }
  virtual ~TrialRigidCluster() {}
  //@}
};

inline std::shared_ptr<TrialRigidCluster> MakeTrialRigidCluster(argtype args = argtype()) {
  return std::make_shared<TrialRigidCluster>(args); }

/**
  Rigid translation of a cluster of particles.
 */
class TrialTranslateCluster : public Trial {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - Trial arguments.
    - SelectCluster arguments.
    - Tunable arguments.
   */
  explicit TrialTranslateCluster(argtype args = argtype());
  explicit TrialTranslateCluster(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialTranslateCluster>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialTranslateCluster>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialTranslateCluster(std::istream& istr);
  virtual ~TrialTranslateCluster() {}
  //@}
};

inline std::shared_ptr<TrialTranslateCluster> MakeTrialTranslateCluster(argtype args = argtype()) {
  return std::make_shared<TrialTranslateCluster>(args); }

/**
  Rigid rotation of a cluster of particles.
 */
class TrialRotateCluster : public Trial {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - Trial arguments.
    - SelectCluster arguments.
    - Tunable arguments.
   */
  explicit TrialRotateCluster(argtype args = argtype());
  explicit TrialRotateCluster(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRotateCluster>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRotateCluster>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRotateCluster(std::istream& istr);
  virtual ~TrialRotateCluster() {}
  //@}
};

inline std::shared_ptr<TrialRotateCluster> MakeTrialRotateCluster(argtype args = argtype()) {
  return std::make_shared<TrialRotateCluster>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_
