
#ifndef FEASST_FFTW_SCATTERING_FFTW_H_
#define FEASST_FFTW_SCATTERING_FFTW_H_

#include <vector>
#include "math/include/position.h"
#include "monte_carlo/include/analyze.h"
#include <complex>
#include <fftw3.h>

namespace feasst {

/**
  Compute the scattering intensity using FFTW (fftw.org) by gridding the system.
  Each site fills a grid of size bin_spacing, but not neighbors, regardless of size.
 */
class ScatteringFFTW : public Analyze {
 public:
  //@{
  /** @name Arguments
    - bin_spacing: maximum bin spacing in each dimension (default: 0.1).
    - bin_per_side: if != -1, ignore bin spacing. (default: -1).
      Instead, set the number of bins per side.
      Powers of two are faster for FFTW (assumes cubic Domain).
    - delta_rho: determines spacing of q values in 3D->1D integeration (default: 1).
  */
  explicit ScatteringFFTW(argtype args = argtype());
  explicit ScatteringFFTW(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("ScatteringFFTW"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<ScatteringFFTW>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<ScatteringFFTW>(args); }
  explicit ScatteringFFTW(std::istream& istr);
  ~ScatteringFFTW();

  //@}
 private:
  bool fftw_initialized_ = false;
  double bin_spacing_;
  int bins_per_side_;
  double delta_rho_;
  int updates_;
  std::vector<int> num_bin_;  // number of bins in each dimension
  std::vector<double> bin_dist_;  // distance of bin in each dimension
  std::vector<double> lower_bound_; // lowest distance bin (half box)
  std::vector<double> fksq_, sqsq_;

  // fftw - not sure if serializable
  double *in_;  // flip
  fftw_complex *out_;
  fftw_plan plan_;

  // return the number of wave vectors (accounting for symmetry along z)
  int num_q_() const {
    return num_bin_[0]*num_bin_[1]*(num_bin_[2]/2+1);
  }

  // resize fftw according to num_bin_
  void resize_fftw_variables_();

  /// fill the in_ grid with points.
  void fill_grid_(const Configuration& config,
    /// If true, only use particle centers (for sq)
    /// Else, account for solid particles.
    const bool centers);

//  /// fill the in_ grid with points from a single site.
//  void fill_grid_site_(
//    const std::vector<int> first_box,
//    const double diameter,
//    const Position& center);
};

inline std::shared_ptr<ScatteringFFTW> MakeScatteringFFTW(
    argtype args = argtype()) {
  return std::make_shared<ScatteringFFTW>(args);
}

}  // namespace feasst

#endif  // FEASST_FFTW_SCATTERING_FFTW_H_
