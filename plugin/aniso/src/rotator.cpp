
#include <chrono> // sleep
#include <thread> // sleep
#include "utils/include/serialize.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/formula.h"
#include "math/include/golden_search.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/rotator.h"

namespace feasst {

Rotator::Rotator(argtype * args) {
  tmp_3vec_.resize(3);
  unique_tolerance_ = dble("unique_tolerance", args, 1e-5);
  contact_tolerance_ = dble("contact_tolerance", args, 1e-4);
  num_proc_ = integer("num_proc", args, 1);
  proc_ = integer("proc", args, 0);
}
Rotator::Rotator(argtype args) : Rotator(&args) {
  FEASST_CHECK_ALL_USED(args);
}

bool Rotator::ior_in_proc(const int ior) const {
  if (ior % num_proc_ == proc_) {
    return true;
  } else {
    return false;
  }
}

void Rotator::gen_orientations(const int num_orientations_per_pi, const Configuration& config) {
  eulers_.clear();
  stheta_.clear();
  sphi_.clear();
  const double dt = PI/static_cast<double>(num_orientations_per_pi);
  double s1min = 0.; // if domain1==domain2, avoid x < 0 for i/j swap symmetry;
  if (config.type_to_file_name(0) != config.type_to_file_name(1)) {
      s1min = -PI;
  }
  int ior = 0;
  for (double s1 = s1min; s1 < PI + dt/2; s1 += dt) {
    for (double s2 = 0.; s2 < PI + dt/2; s2 += dt) { // phi
      for (double e1 = -PI; e1 < PI + dt/2; e1 += dt) {
        for (double e2 = 0; e2 < PI + dt/2; e2 += dt) {
          for (double e3 = -PI; e3 < PI + dt/2; e3 += dt) {
            if (ior_in_proc(ior)) {
              Euler euler(e1, e2, e3);
              eulers_.push_back(euler);
              stheta_.push_back(s1);
              sphi_.push_back(s2);
            }
            ++ior;
          }
        }
      }
    }
  }
  num_orientations_all_proc_ = ior;
  last_three_.clear();
  last_three_.resize(3);
  for (Position& pos : last_three_) {
    pos.set_to_origin(3);
  }
  last_three_sites_.clear();
  last_three_sites_.resize(eulers_.size());
  for (std::vector<Position>& three : last_three_sites_) {
    three.resize(3);
    for (Position& pos : three) {
      pos.set_to_origin(3);
    }
  }
  unique_.clear();
  unique_.resize(eulers_.size());
  contact_.clear();
  contact_.resize(eulers_.size());
  tmp1_.set_to_origin(3);
  tmp2_.set_to_origin(3);
}

void Rotator::init(System * system, const std::string xyz_file,
    const std::string contact_xyz_file) {
  select_ = MakeTrialSelectParticle({{"particle_type", "1"}});
  select_->select_particle(1, system->configuration());
  select_->set_mobile_original(system);
  rotate_ = MakePerturbRotate();
  translate_ = MakePerturbTranslate();
  origin_ = MakePosition();
  origin_->set_to_origin(3);
  rot_mat_.set_size(3, 3);
  xyz_ = FileXYZ({{"group", "mobile"}, {"append", "true"}});
  contact_f_ = FileXYZ({{"group", "mobile"}, {"append", "true"}});
  FileXYZ xyz_fixed_({{"group", "fixed"}, {"append", "false"}});
  if (!xyz_file.empty()) {
    xyz_file_name_ = xyz_file+".xyz";
    std::ofstream ofs;
    ofs.open(xyz_file_name_, std::ofstream::out | std::ofstream::trunc);
    ofs.close();
  }
  if (!contact_xyz_file.empty()) {
    contact_xyz_file_name_ = contact_xyz_file+".xyz";
    std::ofstream ofs;
    ofs.open(contact_xyz_file_name_, std::ofstream::out | std::ofstream::trunc);
    ofs.close();
  }
  if (!xyz_file.empty() || !contact_xyz_file.empty()) {
    //DEBUG("first site " << system->configuration().particle(0).site(0).position().str());
    xyz_fixed_.write(xyz_file+"_fixed.xyz", system->configuration());
  }
}

void Rotator::set_last_three_sites(const int ior, System * system) {
  const Particle& part = system->configuration().particle(1);
  const int num_sites = part.num_sites();
  for (int site = 0; site < 3; ++site) {
    const Position& site_pos = part.site(num_sites - 1 - site).position();
    for (int dim = 0; dim < 3; ++dim) {
      last_three_sites_[ior][site].set_coord(dim, site_pos.coord(dim));
    }
  }
}

void Rotator::check_last_three_sites(const int ior, System * system) {
  tmp1_.set_from_cartesian(last_three_sites_[ior][1].coord());
  tmp2_.set_from_cartesian(last_three_sites_[ior][2].coord());
  tmp1_.subtract(last_three_sites_[ior][0]);
  tmp2_.subtract(last_three_sites_[ior][0]);
  const double cos = tmp1_.cosine(tmp2_);
  ASSERT(std::abs(cos + 1) > NEAR_ZERO && std::abs(cos - 1) > NEAR_ZERO,
    "The last three sites in the molecules form a line and are therefore " <<
    "insufficient for determining unique orientations.");
}

void Rotator::determine_if_unique(const int ior, const std::vector<int>& iors, const int num_threads, System * system) {
//  unique_[ior] = -1;
  int unique = -1;
  //double x = 0.;
  DEBUG("ior " << ior << " unique " << unique);
  for (int past_ior = 0; past_ior < ior; ++past_ior) {
    DEBUG("past_ior " << past_ior);
    const int t = past_ior % num_threads;
    DEBUG("t " << t);
    //while (past_ior > iors[t]) {
    #ifdef _OPENMP
      while (unique_[past_ior] == -2) {
        //x += unique_[past_ior];
        //std::this_thread::sleep_for(std::chrono::milliseconds(1));
        omp_get_thread_num();
        //INFO("waiting ior " << ior);
      }
    #endif // _OPENMP
    DEBUG("unique_[past_ior] " << unique_[past_ior]);
    if (unique_[past_ior] == -1) {
      DEBUG("last three past " << last_three_sites_[past_ior][0].str() << " "
                              << last_three_sites_[past_ior][1].str() << " "
                              << last_three_sites_[past_ior][2].str());
      DEBUG("last three " << last_three_sites_[ior][0].str() << " "
                         << last_three_sites_[ior][1].str() << " "
                         << last_three_sites_[ior][2].str());
      DEBUG("is eq " <<
        last_three_sites_[past_ior][0].is_equal(last_three_sites_[ior][0], unique_tolerance_) << " " <<
        last_three_sites_[past_ior][1].is_equal(last_three_sites_[ior][1], unique_tolerance_) << " " <<
        last_three_sites_[past_ior][2].is_equal(last_three_sites_[ior][2], unique_tolerance_)
      );
      if (last_three_sites_[past_ior][0].is_equal(last_three_sites_[ior][0], unique_tolerance_) &&
          last_three_sites_[past_ior][1].is_equal(last_three_sites_[ior][1], unique_tolerance_) &&
          last_three_sites_[past_ior][2].is_equal(last_three_sites_[ior][2], unique_tolerance_)) {
        DEBUG("found");
        unique = past_ior;
        break;
      }
    }
  }
  DEBUG("ior " << ior << " unique " << unique);
  unique_[ior] = unique;
  DEBUG("unique_[" << ior << "] " << unique_[ior]);

  // only write the xyz for unique orientations
  if (unique_[ior] == -1) {
    if (!xyz_file_name_.empty()) {
      xyz_.write(xyz_file_name_, system->configuration());
    }
//  } else {
//    INFO("non unique stheta:" << stheta_[ior] << " sphi:" << sphi_[ior] <<
//      " ephi:" << eulers_[ior].phi() << " etheta:" << eulers_[ior].theta() << " epsi:" << eulers_[ior].psi());
  }
}

void Rotator::update_xyz(const int ior, const double displacement, System * system) {
  tmp_3vec_[0] = displacement;
  tmp_3vec_[1] = stheta_[ior];
  tmp_3vec_[2] = sphi_[ior];
  DEBUG("sph " << feasst_str(tmp_3vec_));
  if (tmp_3vec_[0] < 1e-8) {
    com1_.set_to_origin(3);
  } else {
    com1_.set_from_spherical(tmp_3vec_);
  }
  DEBUG("com1 " << com1_.str());
  select_->select_particle(1, system->configuration());
  DEBUG("euler " << eulers_[ior].str());
  eulers_[ior].compute_rotation_matrix(&rot_mat_);
  rotate_->set_revert_possible(true, select_.get());
  rotate_->move(*origin_, rot_mat_, system, select_.get());
  translate_->move(com1_, system, select_.get());
  DEBUG("first atom " << system->configuration().particle(1).site(0).position().str());
  DEBUG("first atom fixed " << system->configuration().particle(0).site(0).position().str());
}

void Rotator::revert(System * system) {
  DEBUG("position of first site of mobile before revert " << system->configuration().particle(1).site(0).position().str());
  system->revert(select_->mobile());
  rotate_->revert(system);
  DEBUG("position of first site of mobile after revert " << system->configuration().particle(1).site(0).position().str());
}

int Rotator::num_unique() const {
  int num = 0;
  for (int uniq : unique_) {
    if (uniq == -1) {
      ++num;
    }
  }
  return num;
}

double Rotator::energy(const int ior, const double displacement, System * system, const int ref_potential) {
  update_xyz(ior, displacement, system);
  double en;
//  if (ref_potential == -1) {
//    DEBUG("here");
//  INFO(select_->mobile().num_particles());
  ASSERT(select_->mobile().num_particles() == 1, "err");
  en = system->perturbed_energy(select_->mobile());
//  } else {
//    DEBUG("here");
//    //DEBUG("selp " << select_->mobile().particle_index(0));
//    en = system->reference_energy(select_->mobile(), ref_potential);
//  }
  revert(system);
  DEBUG("position of first site of mobile after revert " << system->configuration().particle(1).site(0).position().str());
  return en;
}

double Rotator::contact_distance(const int ior, System * system, const int ref_potential) {
  if (unique_[ior] != -1) {
    contact_[ior] = contact_[unique_[ior]];
    return contact_[unique_[ior]];
  }
  GoldenSearch minimize({{"tolerance", str(contact_tolerance_)},
    {"lower", "0"},
    {"upper", str(system->configuration().domain().max_side_length()/2)}});
  ContactObjective objective(this, system, ior, ref_potential);
  const double dist = minimize.minimum(&objective) + 2.*contact_tolerance_;
  contact_[ior] = dist;
  if (!contact_xyz_file_name_.empty()) {
    update_xyz(ior, dist, system);
    contact_f_.write(contact_xyz_file_name_, system->configuration());
    revert(system);
  }
  return dist;
}

ContactObjective::ContactObjective(Rotator * rotator, System * system, const double ior, const int ref_pot) {
  rotator_ = rotator;
  system_ = system;
  ior_ = ior;
  ref_pot_ = ref_pot;
}

double ContactObjective::evaluate(const double distance) const {
  const double en = rotator_->energy(ior_, distance, system_, ref_pot_);
  //DEBUG("dist " << distance << " en " << en);
  if (en > 1) {
    return 1e15/(distance + 0.1);
  } else {
    return distance + 0.1;
  }
}

}  // namespace feasst
