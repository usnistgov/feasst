#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/select.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/particle.h"
#include "configuration/include/domain.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"
#include "charge/include/slab_correction.h"

namespace feasst {

FEASST_MAPPER(SlabCorrection, argtype({{"dimension", "0"}}));

SlabCorrection::SlabCorrection(argtype * args) {
  class_name_ = "SlabCorrection";
  dimension_ = integer("dimension", args);
  data_.get_dble_1D()->resize(1);
}
SlabCorrection::SlabCorrection(argtype args) : SlabCorrection(&args) {
  feasst_check_all_used(args);
}

void SlabCorrection::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(2096, ostr);
  //feasst_serialize(conversion_factor_, ostr);
}

SlabCorrection::SlabCorrection(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2096, "unrecognized verison: " << version);
  //feasst_deserialize(&conversion_factor_, istr);
}

class SumDipole : public LoopConfigOneBody {
 public:
  explicit SumDipole(const int dim, const int charge_index)
    : dim_(dim), charge_index_(charge_index) {}
  void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) override {
    const int type = site.type();
    dipole_ += config.model_params().select(charge_index_).value(type)*
               site.position().coord(dim_);
  }
  double dipole() const { return dipole_; }
 private:
  int dim_;
  int charge_index_;
  double dipole_ = 0.;
};

double SlabCorrection::net_dipole(const Configuration& config) const {
  SumDipole sum(dimension_, charge_index());
  VisitConfiguration().loop(config, &sum);
  return sum.dipole();
}

//class SumCharge : public LoopConfigOneBody {
// public:
//  explicit SumCharge() {}
//  void work(const Site& site,
//      const Configuration& config,
//      const LoopDescriptor& data) override {
//    const int type = site.type();
//    charge_ += config.model_params().charge().value(type);
//  }
//  double charge() const { return charge_; }
// private:
//  double charge_ = 0.;
//};

//double SlabCorrection::net_charge(const Configuration& config,
//    const Select& select) const {
//  double charge = 0;
//  for (int sp = 0; sp < select.num_particles(); ++sp) {
//    for (int ss : select.site_indices(sp)) {
//      const int type = config.select_particle(select.particle_index(sp)).site(ss).type();
//      charge += config.model_params().charge().value(type);
//    }
//  }
//  return charge;
//}

double SlabCorrection::dipole_to_en(const double dipole,
    const Configuration& config,
    const ModelParams& params) const {
  const double conv = params.constants().charge_conversion();
  const double volume = config.domain().volume();
  return conv*2*PI/volume*dipole*dipole;
}

void SlabCorrection::compute(ModelOneBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  dipole_new_ = net_dipole(*config);
  DEBUG("dipole_new " << dipole_new_);
  stored_energy_new_ = dipole_to_en(dipole_new_, *config, model_params);
  DEBUG("stored_energy_ " << stored_energy_new_);
  set_energy(stored_energy_new_);
  finalizable_ = true;
}

void SlabCorrection::compute(ModelOneBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  // compute the dipole of the selection
  double sel_dipole = 0;
  for (int sel_part = 0; sel_part < selection.num_particles(); ++sel_part) {
    const int part_index = selection.particle_index(sel_part);
    const Particle& part = config->select_particle(part_index);
    for (int site_index : selection.site_indices(sel_part)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        sel_dipole += model_params.select(charge_index()).value(site.type())*
                      site.position().coord(dimension_);
      }
    }
  }

  // compute new dipole by adding or subtracting dipole of selection
  ASSERT(group_index == 0, "not implemented");
  const int state = selection.trial_state();
  DEBUG("state " << state);
  DEBUG("sel_dipole " << sel_dipole);
  if (state == 0) {
    dipole_new_ = dipole_ - sel_dipole;
  } else if (state == 1) {
    dipole_new_ += sel_dipole;
  } else if (state == 2) {
    dipole_new_ = dipole_ - sel_dipole;
  } else if (state == 3) {
    dipole_new_ = dipole_ + sel_dipole;
  } else {
    FATAL("unrecognized state: " << state);
  }

  // compute new energy
  if (state != 0) {
    stored_energy_new_ = dipole_to_en(dipole_new_, *config, model_params);
  }

  // compute energy change
  double enrg;
  if (state == 0) {
    enrg = stored_energy();
  } else if (state == 1) {
    enrg = stored_energy_new_;
  } else if (state == 2) {
    enrg = stored_energy() - stored_energy_new_;
  } else if (state == 3) {
    enrg = stored_energy_new_ - stored_energy();
  } else {
    FATAL("unrecognized state: " << state);
  }
  DEBUG("enrg: " << enrg);
  DEBUG("stored_energy_ " << stored_energy() << " "
       "stored_energy_new_ " << stored_energy_new_);
  set_energy(enrg);
  finalizable_ = true;
}

void SlabCorrection::finalize(const Select& select, Configuration * config) {
  VisitModel::finalize(select, config);
  dipole_ = dipole_new_;
  *stored_energy_() = stored_energy_new_;
  finalizable_ = false;
}

}  // namespace feasst
