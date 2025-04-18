
#ifndef FEASST_STEPPERS_PAIR_DISTRIBUTION_H_
#define FEASST_STEPPERS_PAIR_DISTRIBUTION_H_

#include "math/include/histogram.h"
#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"
#include "system/include/visit_model_intra.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Compute PairDistribution from the inner two-body vistor
 */
class PairDistributionInner : public ModelTwoBody {
 public:
  PairDistributionInner();

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::vector<std::vector<Histogram> > radial_;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<PairDistributionInner>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit PairDistributionInner(std::istream& istr);
  virtual ~PairDistributionInner();

 protected:
  void serialize_pair_distribution_inner_(std::ostream& ostr) const;
};

typedef std::pair<double, std::vector<std::vector<double> > > grbintype;
typedef std::vector<grbintype> grtype;

/**
  Compute pair distributions.
  This class effectively functions as an Analyze because it does not change the
  System.
  However, the pair distance calculation utilizes Potential::energy, which
  requires a pointer to Configuration.
 */
class PairDistribution : public Modify {
 public:
  //@{
  /** @name Arguments
    - dr: radial distribution bin size (default: 0.1).
    - print_intra: print the intramolecular distributions (default: false).
   */
  explicit PairDistribution(argtype args = argtype());
  explicit PairDistribution(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(MonteCarlo * mc) override;
  std::string write(MonteCarlo * mc) override;

  const grtype& radial(const Configuration& config);

  //const std::vector<std::vector<Histogram> >& radial() const { return radial_; }

  // serialize
  std::string class_name() const override { return std::string("PairDistribution"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<PairDistribution>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<PairDistribution>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit PairDistribution(std::istream& istr);
  explicit PairDistribution(const Modify& pair_distribution);

  //@}
 private:
  double dr_;
  bool print_intra_;
  VisitModel inter_visit_;
  PairDistributionInner inter_;
  VisitModelIntra intra_visit_;
  PairDistributionInner intra_;
  ModelParams params_;
  int num_updates_;

  // temporary and not serialized
  grtype radial_;
};

inline std::shared_ptr<PairDistribution> MakePairDistribution(
    argtype args = argtype()) {
  return std::make_shared<PairDistribution>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_PAIR_DISTRIBUTION_H_
