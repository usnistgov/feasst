
#ifndef FEASST_EWALD_EWALD_H_
#define FEASST_EWALD_EWALD_H_

#include <vector>
#include <sstream>
#include <algorithm>
#include "configuration/include/configuration.h"
#include "system/include/visit_model.h"

namespace feasst {

// HWH Update: use a finalize-heavy instead of revert-heavy strategy
// This means that new eik are calculated and stored separate from system wide properties,
// and system-wide only udpated upon finalize
/**
  This implementation appears to work correctly.
  However, it is not optimized and will be subject to change in future versions.

  HWH: Add Ewald citations.

  Ewald is not currently supported for use as a reference state for dual-cut
  configurational bias.
  While Ewald is supported for use as the full potential for dual-cut
  configurational bias, Ewald cannot and should not be used as one of the
  reference states because (1) it is slow and (2) see Ewald::compute() with
  selections.

  Note to HWH:
  For reverting, there are a few aspects
  1. The eik per-site properties are not reverted by the selection, so must be reverted manually
  2. The structure factor must also be reverted, via perturb->system->potential interface.
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
  // HWH move this to private? and put within wave vectory storage? (or separate from box changes?)
  void update_wave_vectors(const Configuration& config);

  /// Initialize custom site properties to store wave vector components
  void init_wave_vector_storage(Configuration * config, const int group_index = 0);
  void init_wave_vector_storage(Configuration * config, const Select& selection);

  /// Update the site-stored eik properties of the selection, and also the
  /// structure factor.
  void update_struct_fact_eik(const Select& selection,
                              Configuration * config,
                              std::vector<double> * struct_fact_real,
                              std::vector<double> * struct_fact_imag);

  /// Add "eik" to list of config's excluded properties during selection updates
  // HWH do the same for VisitModelCell
  // HWH add init_wave_vector_storage here
  void precompute(Configuration * config) override {
    config->add_excluded_property("eik");
  }

  /// Compute interactions of entire group in configuration from scratch.
  /// This is not optimized for smaller perturbations to the configuration.
  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override {
    // for entire configuration, set stored previous energy to zero
    // store_energy_struct_fact_();
    std::fill(struct_fact_real_new_.begin(), struct_fact_real_new_.end(), 0.);
    std::fill(struct_fact_imag_new_.begin(), struct_fact_imag_new_.end(), 0.);
    update_struct_fact_eik(config->group_select(group_index), config,
                           &struct_fact_real_new_,
                           &struct_fact_imag_new_);
    const double conversion = model_params.constants()->charge_conversion();
    stored_energy_new_ = conversion*fourier_energy_(struct_fact_real_new_,
                                                    struct_fact_imag_new_);
    DEBUG("stored_energy_ " << stored_energy_new_);
    set_energy(stored_energy_new_);
  }

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
    This assumes that there will be a follow up of exactly one 1, "move"
    state which will update the eik.
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
      const ModelOneBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override {
    ASSERT(group_index == 0, "group index cannot be varied because redundant." <<
      "otherwise implement filtering of selection based on group.");
    double enrg = 0.;
    DEBUG("selection.trial_state() " << selection.trial_state());
    ASSERT(selection.trial_state() == 0 ||
           selection.trial_state() == 1 ||
           selection.trial_state() == 2 ||
           selection.trial_state() == 3,
      "unrecognized trial_state: " << selection.trial_state());

    // initialize new structure factor, unless its a new move position
    if (selection.trial_state() != 1) {
      struct_fact_real_new_ = struct_fact_real_;
      DEBUG("size " << struct_fact_real_new_.size() << " " << struct_fact_real_.size());
      struct_fact_imag_new_ = struct_fact_imag_;
    }

    // if "old" half of move, store eik for reverting
    if (selection.trial_state() == 0) {
      // check and resize
      if (static_cast<int>(old_eiks_.size()) != selection.num_particles()) {
        old_eiks_.resize(selection.num_particles());
      }
      for (int ipart = 0; ipart < selection.num_particles(); ++ipart) {
        std::vector<Properties> * eiks = &old_eiks_[ipart];
        if (static_cast<int>(eiks->size()) != selection.num_sites(ipart)) {
          eiks->resize(selection.num_sites(ipart));
        }
        for (int isite = 0; isite < selection.num_sites(ipart); ++isite) {
          Properties * eik = &(*eiks)[isite];
          const int part_index = selection.particle_index(ipart);
          const int site_index = selection.site_index(ipart, isite);
          (*eik) = config->select_particle(part_index).site(site_index).properties();
        }
      }
      revertable_ = true;
      old_config_ = config;
      DEBUG("setting revertable");
    } else if (selection.trial_state() != 1) {
      revertable_ = false;
    }

    update_struct_fact_eik(selection, config, &struct_fact_real_new_,
                                              &struct_fact_imag_new_);
    // compute new energy
    if (selection.trial_state() != 0) {
      const double conversion = model_params.constants()->charge_conversion();
      stored_energy_new_ = conversion*fourier_energy_(struct_fact_real_new_,
                                                      struct_fact_imag_new_);
    }
    if (selection.trial_state() == 0) {
      enrg = stored_energy_;
    } else if (selection.trial_state() == 1) {
      enrg = stored_energy_new_;
    } else if (selection.trial_state() == 2) {
      enrg = stored_energy_ - stored_energy_new_;
    } else if (selection.trial_state() == 3) {
      enrg = stored_energy_new_ - stored_energy_;
    }
    DEBUG("enrg: " << enrg);
    DEBUG("stored_energy_ " << stored_energy_ << " "
          "stored_energy_new_ " << stored_energy_new_);
    set_energy(enrg);
  }

  void revert(const Select& select) override {
    if (revertable_) {
      DEBUG("reverting");
      for (int ipart = 0; ipart < select.num_particles(); ++ipart) {
        const int part_index = select.particle_index(ipart);
        for (int isite = 0; isite < select.num_sites(ipart); ++isite) {
          const int site_index = select.site_index(ipart, isite);
          const Site& site = old_config_->select_particle(part_index).site(site_index);
          const int eikrx0_index = find_eikrx0_(site);
          const std::vector<double>& vals = old_eiks_[ipart][isite].values();
          for (int iprop = 0; iprop < static_cast<int>(vals.size()); ++iprop) {
            old_config_->set_site_property(eikrx0_index + iprop, vals[iprop], part_index, site_index);
          }
        }
      }
    }
//    ERROR("shouldn't be here");
//    struct_fact_real_ = struct_fact_real_old_;
//    struct_fact_imag_ = struct_fact_imag_old_;
//    DEBUG("reverting, stored_energy_ " << stored_energy_);
//    stored_energy_ = stored_energy_old_;
  }

  // HWH refactor Ewald for finalization (e.g., do not enter eiks until finalize?)
  void finalize(const Select& select) override {
    DEBUG("finalizing");
    ASSERT(struct_fact_real_new_.size() > 0, "hi");
    stored_energy_ = stored_energy_new_;
    struct_fact_real_ = struct_fact_real_new_;
    struct_fact_imag_ = struct_fact_imag_new_;
    revertable_ = false;
  }

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
  const int dimension_ = 3;
  double stored_energy_ = 0.;

  std::vector<double> struct_fact_real_;
  std::vector<double> struct_fact_imag_;

  // temporary
  std::vector<double> struct_fact_real_new_;
  std::vector<double> struct_fact_imag_new_;
  double stored_energy_new_ = 0.;

  // temporary
  // HWH not sure this is the best way to store and revert eiks
  // HWH but a refactor would require argument SelectParticle * select
  // HWH to put eiks in selection and no longer exclude them from update.
  std::vector<std::vector<Properties> > old_eiks_;
  bool revertable_ = false;
  Configuration * old_config_;

//  void store_energy_struct_fact_() {
//    stored_energy_old_ = stored_energy_;
//    struct_fact_real_old_ = struct_fact_real_;
//    struct_fact_imag_old_ = struct_fact_imag_;
//  }

  std::vector<std::string> eik_gen_();

  double fourier_energy_(const std::vector<double>& struct_fact_real,
                         const std::vector<double>& struct_fact_imag) {
    double en = 0;
    for (int k = 0; k < num_vectors(); ++k) {
      en += wave_prefactor_[k]*(struct_fact_real[k]*struct_fact_real[k]
                              + struct_fact_imag[k]*struct_fact_imag[k]);
    }
    return en;
  }

  double sign_(const Select& select) {
    if (select.trial_state() == 0 || select.trial_state() == 2) {
      return -1.0;
    }
    return 1.0;
  }

  // temporary
  std::string eikrx0_str_ = "eikrx0";
  int find_eikrx0_(const Site& site) {
    int eikrx0_index = 0;
    ASSERT(
      find_in_list(eikrx0_str_, site.properties().names(), &eikrx0_index),
      "eikrx0 doesn't exist");
    return eikrx0_index;
  }
};

inline std::shared_ptr<Ewald> MakeEwald() {
  return std::make_shared<Ewald>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_EWALD_H_
