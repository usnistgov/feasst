
#ifndef FEASST_EWALD_EWALD_H_
#define FEASST_EWALD_EWALD_H_

#include <vector>
#include <sstream>
#include <algorithm>
#include "configuration/include/configuration.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  The Ewald summation accounts for the long-range nature of the electrostatic
  interaction by applying a Gaussian screening charge, computing a Fourier-space
  long-range component, and then correcting for the various spurious terms
  that are included in the Fourier summation, such as self and intra-particle.

  See "Computer simulation of liquids" by M. P. Allen and D. J. Tildesley.
  The LAMMPS and DL_POLY manuals also include thorough descriptions of the Ewald
  summation.

  Ewald is not supported for use as a reference state for dual-cut.
 */
class Ewald : public VisitModel {
 public:
  /**
    args:
    - tolerance: determine the alpha parameter and number of wave vectors by
      specifying the accuracy relative to the energy of two unit charges
      separated by a distance of one unit.
    - tolerance_num_sites: for setting parameters with the tolerance,
      optionally set the number of sites to be used for the parameter
      calculation rather than the currently existing number of sites (default).
    - alpha: optionally specify the alpha parameter in units of inverse length.
    - kxmax: optionally specify the maximum wave vectors in the first dimension.
    - kymax: same as above, but in the second dimension.
    - kzmax: same as above, but in the third dimension.
    - kmax_squared: optionally set the squared maximum integer wave vector for
      cubic domains only, which also sets kxmax, etc.
   */
  Ewald(argtype args = argtype());

  /**
    Recommend an alpha parameter for Ewald as described and implemented in LAMMPS
    https://lammps.sandia.gov/doc/kspace_style.html
    https://doi.org/10.1080/08927029208049126
    https://doi.org/10.1063/1.470043
   */
  void tolerance_to_alpha_ks(
    /// Estimated tolerance for the energy.
    /// This is relative to the energy exerted by two unit charges
    /// separated by a distance of one unit.
    const double tolerance,
    const Configuration& config,
    /// Return the alpha
    double * alpha,
    /// Return the maximum number of Fourier vectors in the x dimension.
    int * kxmax,
    /// as above but in y
    int * kymax,
    /// as above but in z
    int * kzmax);

  // HWH move this to private? (or separate from box changes?)
  /// Precompute the wave vectors within cutoff, coefficients, and also resize
  /// the structure factors.
  void update_wave_vectors(const Configuration& config);

  /// Compute new eiks and update the given structure factor.
  void update_struct_fact_eik(const Select& selection,
    const Configuration& config,
    std::vector<double> * struct_fact_real,
    std::vector<double> * struct_fact_imag,
    std::vector<std::vector<std::vector<double> > > * eik_new) const;

  /// Process tolerance arguments and initialize wave vectors.
  void precompute(Configuration * config) override;

  /// Compute interactions of entire group in configuration from scratch.
  /// This is not optimized for smaller perturbations to the configuration.
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override;

  /**
    Compute by selection is optimized for perturbations to the system.
    Contrary to short range interactions, in this case the energy of the
    entire system is computed from the structure factor.
    Thus, to get an energy contribution of a selection, one must take a
    difference between the configuration with and without the selection's
    contributions to the structure factor.

    For MonteCarlo trials, this means that the type of trial must be taken
    into consideration, which is provided via Select::trial_state().

    For 0, "old", first prepare to return the (previously computed) energy of
    the entire system, and then remove the structure factor contributions of
    the selection, without updating the selection's eik.
    This assumes that there will be a follow up of exactly one state 1, "move"
    which will update the eik.
    Ewald is not currently supported with more than one step in TrialStage.

    For 1, "move", eik are updated and their contributions are added to the
    structure factor.
    Return the new system energy.

    For 2, "remove", remove the eik contributions to the structure factor
    without updating the eik.
    Return the old minus new energy.

    For 3, "add", eik are updated and their contributions are added to the
    structure factor.
    Return the new minus old energy.
   */
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;

  void change_volume(const double delta_volume, const int dimension) override;

  // update structure factors and eiks based on new calculations.
  void finalize(const Select& select, Configuration * config) override;

  int num_vectors() const { return static_cast<int>(wave_prefactor_.size()); }
  int num_kx() const { return num_kx_; }
  int num_ky() const { return num_ky_; }
  int num_kz() const { return num_kz_; }
  int kxmax() const { return kxmax_; }
  int kymax() const { return kymax_; }
  int kzmax() const { return kzmax_; }
  void check_size() const;
  const std::vector<double>& struct_fact_real() const {
    return data_.dble_2D()[0]; }
  const std::vector<double>& struct_fact_imag() const {
    return data_.dble_2D()[1]; }

  /// Return the net charge of the configuration.
  double net_charge(const Configuration& config) const;

  // Return the eik vectors directly.
  const std::vector<std::vector<std::vector<double> > >& eik() const {
    //return const_cast<std::vector<std::vector<std::vector<double> > >&>(*eik_()); }
    return manual_data_.dble_3D(); }

  void check(const Configuration& config) const override;

  void synchronize_(const VisitModel& visit, const Select& perturbed) override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  Ewald(std::istream& istr);
  void serialize(std::ostream& ostr) const override;

 private:
  // HWH serialize
  std::shared_ptr<double> tolerance_, alpha_arg_;
  std::shared_ptr<int> tolerance_num_sites_, kxmax_arg_, kymax_arg_, kzmax_arg_, kmax_sq_arg_;
  int kxmax_, kymax_, kzmax_;
  double kmax_squared_;
  int num_kx_;
  int num_ky_;
  int num_kz_;
  std::vector<double> wave_prefactor_;
  std::vector<int> wave_num_;
  const int dimension_ = 3;
  //double stored_energy_ = 0.;

  // synchronization data
  double stored_energy() const { return data_.dble_1D()[0]; }
  double * stored_energy_() { return &((*data_.get_dble_1D())[0]); }
  std::vector<double> * struct_fact_real_();
  std::vector<double> * struct_fact_imag_();

  // new eik implementation, Ewald contains all eik information.
  // eik_[particle_index][site_index][eik_index]
  std::vector<std::vector<std::vector<double> > > * eik_();
  // temporary
  std::vector<std::vector<std::vector<double> > > eik_new_;

  // not temporary (for sizing)
  std::vector<double> struct_fact_real_new_;
  std::vector<double> struct_fact_imag_new_;
  // temporary
  double stored_energy_new_ = 0.;

  // temporary
  bool finalizable_ = false;

  std::vector<std::string> eik_gen_();

  /// Return the sum of the squared charge.
  double sum_squared_charge_(const Configuration& config) {
    double sum_sq_q = 0.;
    const std::vector<int> num_sites_of_type = config.num_sites_of_type();
    for (int type = 0;
         type < static_cast<int>(num_sites_of_type.size());
         ++type) {
      const double charge = config.model_params().charge().value(type);
      sum_sq_q += charge*charge*num_sites_of_type[type];
    }
    return sum_sq_q;
  }

  /// Return the Fourier root mean squared accuracy for a given dimension.
  double fourier_rms_(
      const double alpha,
      const int kmax,
      const Configuration& config,
      const int dimen,
      const int num_sites);

  int estimate_kmax_(
      const double alpha,
      const Configuration& config,
      const double tolerance,
      const int dimen,
      const int num_sites);

  double fourier_energy_(const std::vector<double>& struct_fact_real,
                         const std::vector<double>& struct_fact_imag);

  double sign_(const Select& select, const int pindex) const;

  void resize_eik_(const Configuration& config);
};

inline std::shared_ptr<Ewald> MakeEwald(argtype args = argtype()) {
  return std::make_shared<Ewald>(args);
}

}  // namespace feasst

#endif  // FEASST_EWALD_EWALD_H_
