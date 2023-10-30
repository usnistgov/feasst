
#ifndef FEASST_CLUSTER_TRIAL_TRANSLATE_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_TRANSLATE_CLUSTER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Rigid translation of a cluster of particles.
 */
class TrialTranslateCluster : public Trial {
 public:
  //@{
  /** @name Arguments
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

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_TRANSLATE_CLUSTER_H_
