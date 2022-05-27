
#ifndef FEASST_CLUSTER_ANALYZE_CLUSTER_H_
#define FEASST_CLUSTER_ANALYZE_CLUSTER_H_

#include "math/include/accumulator.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Harvest cluster size information from any Trial with SelectCluster in its first stage (e.g., TrialTranslateCluster and TrialRotateCluster).
 */
class AnalyzeCluster : public Analyze {
 public:
  explicit AnalyzeCluster(argtype args = argtype());
  explicit AnalyzeCluster(argtype * args);

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  /// Return the cluster size.
  const Accumulator& cluster_size() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("AnalyzeCluster"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeCluster>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<AnalyzeCluster>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit AnalyzeCluster(std::istream& istr);
  explicit AnalyzeCluster(const Analyze& energy);
};

inline std::shared_ptr<AnalyzeCluster> MakeAnalyzeCluster(argtype args = argtype()) {
  return std::make_shared<AnalyzeCluster>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ANALYZE_CLUSTER_H_
