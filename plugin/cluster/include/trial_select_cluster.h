
#ifndef FEASST_CLUSTER_TRIAL_SELECT_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_SELECT_CLUSTER_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "system/include/cluster_criteria.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/**
  Select a random particle.
  Find all particles which interact with any of the sites on the selection
  via full potential (default) or a given reference potential.
 */
class TrialSelectCluster : public TrialSelect {
 public:
  TrialSelectCluster(
    std::shared_ptr<ClusterCriteria> cluster_criteria,
    const argtype& args = argtype());

  /// Return a cluster as selection given one particle in the system.
  void select_cluster(const int first_particle, const System& system,
                      SelectList * select);

  // Same as above, except put selection in mobile.
  void select_cluster(const int first_particle, const System& system) {
    select_cluster(first_particle, system, &mobile_); }

  /// Return a list of selections representing individual cluster.
  std::vector<SelectList> select_clusters(const System& system);

  bool select(const Select& perturbed,
              System* system,
              Random * random) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectCluster(std::istream& istr);
  virtual ~TrialSelectCluster() {}

 protected:
  void serialize_trial_select_cluster_(std::ostream& ostr) const;

 private:
  std::shared_ptr<ClusterCriteria> cluster_criteria_;
  std::shared_ptr<TrialSelectParticle> select_particle_;
};

inline std::shared_ptr<TrialSelectCluster> MakeTrialSelectCluster(
    std::shared_ptr<ClusterCriteria> cluster_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectCluster>(cluster_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_SELECT_CLUSTER_H_
