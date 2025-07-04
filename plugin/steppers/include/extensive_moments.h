
#ifndef FEASST_STEPPERS_EXTENSIVE_MOMENTS_H_
#define FEASST_STEPPERS_EXTENSIVE_MOMENTS_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate extensive moments for derivatives from fluctuations.
  This is currently implemented for a grand canonical ensemble:

  \f$N_i^j N_k^m U^p\f$

  where \f$i,k\f$ are particle types, and the powers \f$j,m,p\f$ are collected
  from 0 to a maximum order cutoff.
 */
class ExtensiveMoments : public Analyze {
 public:
  //@{
  /** @name Arguments
    - max_order: maximum order cutoff for the powers (default: 3).
    - Stepper arguments.
   */
  explicit ExtensiveMoments(argtype args = argtype());
  explicit ExtensiveMoments(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  /// Return the extensive moments
  const Accumulator& moments(const int p, const int m, const int k,
      const int j, const int i) const {
    return moments_[p][m][k][j][i]; }

  // serialize
  std::string class_name() const override { return std::string("ExtensiveMoments"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<ExtensiveMoments>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<ExtensiveMoments>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ExtensiveMoments(std::istream& istr);
  explicit ExtensiveMoments(const Analyze& extensive_moments);

  //@}
 private:
  int max_order_;
  // moments_[p][m][k][j][i]
  std::vector<std::vector<std::vector<std::vector<std::vector<Accumulator> > > > > moments_;
  std::vector<double> u_p_;
  std::vector<std::vector<double> > n_i_j_;
};

inline std::shared_ptr<ExtensiveMoments> MakeExtensiveMoments(argtype args = argtype()) {
  return std::make_shared<ExtensiveMoments>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_EXTENSIVE_MOMENTS_H_
