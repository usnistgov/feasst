
#ifndef FEASST_EWALD_SLAB_CORRECTION_H_
#define FEASST_EWALD_SLAB_CORRECTION_H_

#include "system/include/visit_model.h"

namespace feasst {

/**
  For a periodically repeated slab 3 to 5 times larger than the thickness
  of the slab,

  \f$U = \frac{2\pi}{V} M_z^2\f$

  \f$M_z = \sum^N_{i=i}q_i z_i\f$

  where z is the position in the slab dimension, q is the charge and V is the
  volume.

  I.C. Yeh and M.L. Berkowitz, J. Chem. Phys., 111, 3155, 1999.

  E.R. Smith, Proc. R. Soc. London A, 375, 475-505, 1981.

  P.S. Crozier, R.L. Rowley, E. Spohr and D. Henderson, J. Chem. Phys., 112, 925309257, 2000.

  The non-neutral terms are not yet implemented:

  V. Ballenegger, A. Arnold and J. J. Cerda, J. Chem. Phys. 131, 094107 (2009).
 */
class SlabCorrection : public VisitModel {
 public:
  /**
    args:
    - dimension: dimension of periodically replicated slab.
   */
  explicit SlabCorrection(argtype args);

  //void precompute(Configuration * config) override;

  /// Return the net dipole of the configuration.
  double net_dipole(const Configuration& config) const;

//  /// Return the net charge of the configuration.
//  double net_charge(const Configuration& config,
//    /// The trial state is used to account for sites selected for deletion.
//    const Select& select) const;

  /// Compute interactions of entire group in configuration from scratch.
  /// This is not optimized for smaller perturbations to the configuration.
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override;

  /**
    Compute by selection is optimized for perturbations to the system.
   */
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;

  void finalize(const Select& select, Configuration * config) override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<SlabCorrection>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit SlabCorrection(std::istream& istr);
  virtual ~SlabCorrection() {}

 private:
  double dimension_;
  double dipole_;

  // synchronization data
  double stored_energy() const { return data_.dble_1D()[0]; }
  double * stored_energy_() { return &((*data_.get_dble_1D())[0]); }

  // temporary
  double dipole_new_;
  bool finalizable_ = false;
  double stored_energy_new_ = 0.;

  double dipole_to_en(const double dipole, const Configuration& config,
    const ModelParams& params) const;

  double sum_charge(const Configuration& config, const int state) const;
};

inline std::shared_ptr<SlabCorrection> MakeSlabCorrection(argtype args) {
  return std::make_shared<SlabCorrection>(args);
}

}  // namespace feasst

#endif  // FEASST_EWALD_SLAB_CORRECTION_H_
