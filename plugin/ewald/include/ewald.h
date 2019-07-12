
#ifndef FEASST_EWALD_EWALD_H_
#define FEASST_EWALD_EWALD_H_

#include <vector>
#include <sstream>
#include <algorithm>
#include "system/include/physical_constants.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Note to HWH:
  For reverting, there are a few aspects
  1. The eik per-site properties are reverted by the selection
  2. The structure factor must also be reverted, through perturb->system->potential interface.
  3. Visitor needs to know if its a 'new' or 'old' configuration ? (what about transfers?)

  HWH: If selection = all, don't do difference
  HWH: If selection not all (how know? config?) then calc energy first before updating struct.
 */
class Ewald : public VisitModel {
 public:
  Ewald() {}

  /// Set all k vectors according to similar size and use symmetry for the x.
  void set_kmax_squared(const double kmax_squared);

  /// Precompute the wave vectors within cutoff, coefficients, and also resize
  /// the structure factors.
  void update_wave_vectors(const Configuration& config);

  /// Initialize custom site properties to store wave vector components
  void init_wave_vector_storage(Configuration * config, const int group_index = 0);
  void init_wave_vector_storage(Configuration * config, const Select& selection);

  /// Compute by selection requires subtracting from the total energy,
  /// which is stored for optimization
  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override {
    ASSERT(group_index == 0, "group index cannot be varied because redundant." <<
      "otherwise implement filtering of selection based on group.");
    store_struct_fact_();
    update_eik(selection, config);
    stored_energy_old_ = stored_energy_;
    stored_energy_ = fourier_energy_();
    const double en = sign_(selection)*(stored_energy_ - stored_energy_old_);
    set_energy(en);
  }

  /// Compute by entire selection does not require subtraction as above.
  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override {
    // for entire configuration, set stored previous energy to zero
    store_struct_fact_();
    std::fill(struct_fact_real_.begin(), struct_fact_real_.end(), 0.);
    std::fill(struct_fact_imag_.begin(), struct_fact_imag_.end(), 0.);
    update_eik(config->group_select(group_index), config);
    stored_energy_old_ = stored_energy_;
    stored_energy_ = fourier_energy_();
    set_energy(stored_energy_);
  }

  void revert() override {
    struct_fact_real_ = struct_fact_real_old_;
    struct_fact_imag_ = struct_fact_imag_old_;
    stored_energy_ = stored_energy_old_;
  }

  /// Update the site-stored eik properties and also the structure factor.
  void update_eik(const Select& selection, Configuration * config);

  int num_vectors() const { return static_cast<int>(wave_prefactor_.size()); }

  void check_size() const {
    ASSERT(wave_prefactor_.size() == wave_num_.size(), "size err");
  }

  std::vector<double> struct_fact_real() const { return struct_fact_real_; }
  std::vector<double> struct_fact_imag() const { return struct_fact_imag_; }

  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  Ewald(std::istream& istr);
  void serialize(std::ostream& ostr) const override;

 private:
  const std::string class_name_ = "Ewald";
  int kmax_;
  int kmax_squared_;
  int kxmax_;
  int kymax_;
  int kzmax_;
  std::vector<double> wave_prefactor_;
  std::vector<int> wave_num_;
  std::vector<double> struct_fact_real_;
  std::vector<double> struct_fact_imag_;
  std::vector<double> struct_fact_real_old_;
  std::vector<double> struct_fact_imag_old_;
  const int dimension_ = 3;
  double stored_energy_ = 0.;
  double stored_energy_old_ = 0.;

  void store_struct_fact_() {
    struct_fact_real_old_ = struct_fact_real_;
    struct_fact_imag_old_ = struct_fact_imag_;
  }

  double fourier_energy_() {
    double en = 0;
    for (int k = 0; k < num_vectors(); ++k) {
      en += wave_prefactor_[k]*(struct_fact_real_[k]*struct_fact_real_[k]
                              + struct_fact_imag_[k]*struct_fact_imag_[k]);
    }
    return charge_conversion*en;
  }

  double sign_(const Select& select) {
    if (select.trial_state() == 0) {
      return -1.0;
    }
    return 1.0;
  }
};

inline std::shared_ptr<Ewald> MakeEwald() {
  return std::make_shared<Ewald>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_EWALD_H_
