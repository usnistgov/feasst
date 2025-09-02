
#ifndef FEASST_EWALD_EWALD_H_
#define FEASST_EWALD_EWALD_H_

#include <vector>
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

  Following the description in the classic DL-POLY user manual (version 1.9),
  if the real space basis vectors are given by \f$\vec{a}, \vec{b}, \vec{c}\f$,
  then reciprocal space basis vectors are

  \f$V = \vec{a} \cdot \vec{b} \times \vec{c}\f$

  \f$\vec{u} = 2\pi\vec{b}\times\vec{c}/V\f$

  \f$\vec{v} = 2\pi\vec{b}\times\vec{c}/V\f$

  \f$\vec{w} = 2\pi\vec{b}\times\vec{c}/V\f$

  See <a href="../tutorial/tutorial_0_spce_ref_config.html">this tutorial</a> for an example.
 */
class Ewald : public VisitModel {
 public:
  //@{
  /** @name Arguments
    - tolerance: determine the alpha parameter and number of wave vectors by
      specifying the accuracy relative to the energy of two unit charges
      separated by a distance of one unit.
      As described and implemented in LAMMPS, see:
      https://docs.lammps.org/kspace_style.html,
      https://doi.org/10.1080/08927029208049126,
      https://doi.org/10.1063/1.470043.
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
  explicit Ewald(argtype args = argtype());
  explicit Ewald(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /**
    Recommend an alpha parameter for Ewald as described and implemented in LAMMPS
    https://docs.lammps.org/kspace_style.html
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
  void update_wave_vectors(const Configuration& config,
    const double kmax_squared,
    std::vector<double> * wave_prefactor,
    std::vector<int> * wave_num,
    double* ux, double* uy, double* uz, double* vy, double* vz, double* wz) const;

  void update_kmax_squared_(const Configuration& config, double * kmax_squared) const;
  void resize_struct_fact_new_(const int num_vectors);

  /// Compute new eiks and update the given structure factor.
  void update_struct_fact_eik(const Select& selection,
    const Configuration& config,
    const std::vector<double>& wave_prefactor,
    const std::vector<int>& wave_num,
    const double ux, const double uy, const double uz,
    const double vy, const double vz, const double wz,
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

    For 4, "volume", use the compute without selection.
   */
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;

  void change_volume(const double delta_volume, const int dimension,
    Configuration * config) override;

  // update structure factors and eiks based on new calculations.
  void finalize(const Select& select, Configuration * config) override;

  /// Return the wave vector number for a given dimension.
  int wave_num(const int vector_index, const int dim) const;

  // Return the prefactor in the fourier-space term for each vector.
  //const std::vector<double>& wave_prefactor() const { return wave_prefactor_; }

  /// Return the number of vectors in the x dimension.
  int num_kx() const { return num_kx_; }

  /// Return the number of vectors in the y dimension.
  int num_ky() const { return num_ky_; }

  /// Return the number of vectors in the z dimension.
  int num_kz() const { return num_kz_; }

  /// Return the cutoff of the wave vector in the x dimension.
  int kxmax() const { return kxmax_; }

  /// Return the cutoff of the wave vector in the y dimension.
  int kymax() const { return kymax_; }

  /// Return the cutoff of the wave vector in the z dimension.
  int kzmax() const { return kzmax_; }

  // See synchronize data:
  // Return the spherical cutoff of the wave vectors.
  // double kmax_squared() const { return kmax_squared_; }

//  /// Return the eik vectors directly.
//  double eik(const int part_index, const int site_index,
//    const int vector_index, const int dim, const bool real = true) const;

  /// Same as above, but for all particles, sites and wave vectors.
  const std::vector<std::vector<std::vector<double> > >& eik() const {
    return manual_data_.dble_3D(); }

  /// Return the real part of the structure factor for a given vector index
  /// corresponding with wave_prefactor and wave_num.
  double struct_fact_real(const int vector_index) const {
    return struct_fact_real()[vector_index]; }

  /// Return the imaginary part of the structure factor for a given vector index
  /// corresponding with wave_prefactor and wave_num.
  double struct_fact_imag(const int vector_index) const {
    return struct_fact_imag()[vector_index]; }

  const std::vector<double>& struct_fact_real() const {
    return data_.dble_2D()[0]; }

  /// Return the imaginary part of the structure factor.
  const std::vector<double>& struct_fact_imag() const {
    return data_.dble_2D()[1]; }

  /// Return the net charge of the configuration.
  double net_charge(const Configuration& config) const;

  void check(const Configuration& config) const override;

  void synchronize_(const VisitModel& visit, const Select& perturbed) override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<Ewald>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<Ewald>(args); }
  Ewald(std::istream& istr);
  void serialize(std::ostream& ostr) const override;

  //@}
 private:
  const int dimension_ = 3;
  // HWH serialize
  std::shared_ptr<double> tolerance_;
  std::shared_ptr<double> alpha_arg_;
  std::shared_ptr<int> tolerance_num_sites_, kxmax_arg_, kymax_arg_, kzmax_arg_, kmax_sq_arg_;
  int kxmax_;
  int kymax_;
  int kzmax_;
  int num_kx_;
  int num_ky_;
  int num_kz_;
  double uy_, uz_, vz_;
  double uy_new_, uz_new_, vz_new_;
  std::vector<int> wave_num_;
  std::vector<int> wave_num_new_;

  // synchronization data
  double stored_energy() const { return data_.dble_1D()[0]; }
  double * stored_energy_() { return &((*data_.get_dble_1D())[0]); }
  double * kmax_squared_() { return &(*data_.get_dble_1D())[1]; }
  const double kmax_squared() const { return data_.dble_1D()[1]; }
  double * kmax_squared_new_() { return &(*data_.get_dble_1D())[2]; }
  const double kmax_squared_new() const { return data_.dble_1D()[2]; }
  double * ux_() { return &(*data_.get_dble_1D())[3]; }
  const double ux() const { return data_.dble_1D()[3]; }
  double * vy_() { return &(*data_.get_dble_1D())[4]; }
  const double vy() const { return data_.dble_1D()[4]; }
  double * wz_() { return &(*data_.get_dble_1D())[5]; }
  const double wz() const { return data_.dble_1D()[5]; }
  double * ux_new_() { return &(*data_.get_dble_1D())[6]; }
  const double ux_new() const { return data_.dble_1D()[6]; }
  double * vy_new_() { return &(*data_.get_dble_1D())[7]; }
  const double vy_new() const { return data_.dble_1D()[7]; }
  double * wz_new_() { return &(*data_.get_dble_1D())[8]; }
  const double wz_new() const { return data_.dble_1D()[8]; }
  std::vector<double> * struct_fact_real_() { return &((*data_.get_dble_2D())[0]); }
  std::vector<double> * struct_fact_imag_() { return &((*data_.get_dble_2D())[1]); }
  std::vector<double> * struct_fact_real_new_() { return &((*data_.get_dble_2D())[2]); }
  const std::vector<double>& struct_fact_real_new() const { return data_.dble_2D()[2]; }
  std::vector<double> * struct_fact_imag_new_() { return &((*data_.get_dble_2D())[3]); }
  const std::vector<double>& struct_fact_imag_new() const { return data_.dble_2D()[3]; }
  std::vector<double> * wave_prefactor_() { return &((*data_.get_dble_2D())[4]); }
  const std::vector<double>& wave_prefactor() const { return data_.dble_2D()[4]; }
  std::vector<double> * wave_prefactor_new_() { return &((*data_.get_dble_2D())[5]); }
  const std::vector<double>& wave_prefactor_new() const { return data_.dble_2D()[5]; }

  // Ewald contains all eik information.
  // eik_[particle_index][site_index][eik_index]
  std::vector<std::vector<std::vector<double> > > * eik_();
  // temporary
  std::vector<std::vector<std::vector<double> > > * eik_new_();
  const std::vector<std::vector<std::vector<double> > >& eik_new() const;

  // temporary
  double stored_energy_new_ = 0.;
  bool finalizable_ = false;

  /// Return the sum of the squared charge.
  /// If Configuraiton has no particles but num_sites != 0, assume first particle type.
  double sum_squared_charge_(const Configuration& config, const int num_sites);

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
                         const std::vector<double>& struct_fact_imag,
                         const std::vector<double>& wave_prefactor);

  double sign_(const Select& select, const int pindex) const;

  void resize_eik_(const Configuration& config);
  void resize_eik_(const std::vector<std::vector<std::vector<double> > >& eik2);
};

inline std::shared_ptr<Ewald> MakeEwald(argtype args = argtype()) {
  return std::make_shared<Ewald>(args);
}

}  // namespace feasst

#endif  // FEASST_EWALD_EWALD_H_
