
#ifndef FEASST_CLUSTER_CALCULATE_CLUSTER_H_
#define FEASST_CLUSTER_CALCULATE_CLUSTER_H_

#include "math/include/accumulator.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Compute cluster sizes.
  This class is in development and may not function properly.
 */
class CalculateCluster : public Modify {
 public:
  explicit CalculateCluster(argtype args = argtype());
  explicit CalculateCluster(argtype * args);

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(MonteCarlo * mc) override;
  std::string write(MonteCarlo * mc) override;

  /// Return the cluster size.
  const Accumulator& cluster_size() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("CalculateCluster"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CalculateCluster>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<CalculateCluster>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit CalculateCluster(std::istream& istr);
  explicit CalculateCluster(const Modify& energy);
};

inline std::shared_ptr<CalculateCluster> MakeCalculateCluster(argtype args = argtype()) {
  return std::make_shared<CalculateCluster>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_CALCULATE_CLUSTER_H_
