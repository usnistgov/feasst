#include <cmath>  // abs
#include <chrono> // sleep
#include <thread> // sleep
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/constants.h"
#include "math/include/formula.h"
#include "math/include/golden_search.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/rotator.h"

namespace feasst {

Rotator::Rotator(argtype * args) {
  tmp_3vec_.resize(3);
  unique_tolerance_ = dble("unique_tolerance", args, 1e-5);
  contact_tolerance_ = dble("contact_tolerance", args, 1e-4);
  num_proc_ = integer("num_proc", args, 1);
  proc_ = integer("proc", args, 0);
  hard_limit_u_ = dble("hard_limit_u", args, hard_u_);
}
Rotator::Rotator(argtype args) : Rotator(&args) {
  feasst_check_all_used(args);
}

bool Rotator::ior_in_proc(const int ior) const {
  if (ior % num_proc_ == proc_) {
    return true;
  } else {
    return false;
  }
}

std::vector<std::vector<double> > Rotator::gen_global_bounds(const Configuration& config) const {
  double s1max = 2*PI;
  // if domain1==domain2, avoid x < 0 for i/j swap symmetry;
  DEBUG(config.type_to_file_name(0));
  if (config.num_particle_types() == 1) {
    s1max = PI;
  } else if (config.num_particle_types() == 2) {
    if (config.type_to_file_name(0) == config.type_to_file_name(1)) {
      s1max = PI;
    }
  } else {
    FATAL("unrecognized number of particle types");
  }
  if (config.particle(1).num_sites() == 1) {
    return {{0, s1max}, {0, PI}};
  }
  return {{0, s1max}, {0, PI}, {-PI, PI}, {0, PI}, {-PI, PI}};
}

void Rotator::gen_orientations(const int num_orientations_per_pi, const Configuration& config, std::vector<std::vector<double> > bounds) {
  eulers_.clear();
  stheta_.clear();
  sphi_.clear();
  if (bounds.size() == 0) {
    bounds = gen_global_bounds(config);
  }
  const double dt = (bounds[1][1]-bounds[1][0])/static_cast<double>(num_orientations_per_pi);
  int ior = 0;
  int is1 = 0;
  if (config.particle(1).num_sites() == 1) {
    sizes_.resize(2);
    for (double s1 = bounds[0][0]; s1 < bounds[0][1] + dt/2; s1 += dt) { // theta
      int is2 = 0;
      for (double s2 = bounds[1][0]; s2 < bounds[1][1] + dt/2; s2 += dt) { // phi
        if (ior_in_proc(ior)) {
          stheta_.push_back(s1);
          sphi_.push_back(s2);
          indices_.push_back({is1, is2});
          //INFO(s1 << " " << s2);
        }
        ++ior;
        ++is2;
      }
      sizes_[1] = is2;
      ++is1;
    }
    sizes_[0] = is1;
  } else {
    sizes_.resize(5);
    for (double s1 = bounds[0][0]; s1 < bounds[0][1] + dt/2; s1 += dt) { // theta
      int is2 = 0;
      for (double s2 = bounds[1][0]; s2 < bounds[1][1] + dt/2; s2 += dt) { // phi
        int ie1 = 0;
        for (double e1 = bounds[2][0]; e1 < bounds[2][1] + dt/2; e1 += dt) { // ephi
          int ie2 = 0;
          for (double e2 = bounds[3][0]; e2 < bounds[3][1] + dt/2; e2 += dt) { // etheta
            int ie3 = 0;
            for (double e3 = bounds[4][0]; e3 < bounds[4][1] + dt/2; e3 += dt) { // epsi
              if (ior_in_proc(ior)) {
                Euler euler(e1, e2, e3);
                eulers_.push_back(euler);
                stheta_.push_back(s1);
                sphi_.push_back(s2);
                indices_.push_back({is1, is2, ie1, ie2, ie3});
                //INFO(s1 << " " << s2 << " " << e1 << " " << e2 << " " << e3);
              }
              ++ior;
              ++ie3;
            }
            sizes_[4] = ie3;
            ++ie2;
          }
          sizes_[3] = ie2;
          ++ie1;
        }
        sizes_[2] = ie1;
        ++is2;
      }
      sizes_[1] = is2;
      ++is1;
    }
    sizes_[0] = is1;
  }
  num_orientations_all_proc_ = ior;
  last_three_.clear();
  last_three_.resize(3);
  for (Position& pos : last_three_) {
    pos.set_to_origin(3);
  }
  last_three_sites_.clear();
  last_three_sites_.resize(num_orientations());
  for (std::vector<Position>& three : last_three_sites_) {
    three.resize(3);
    for (Position& pos : three) {
      pos.set_to_origin(3);
    }
  }
  unique_.clear();
  unique_.resize(num_orientations());
  contact_.clear();
  contact_.resize(num_orientations());
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
  ASSERT(system->configuration().particle(0).site(0).position().squared_distance(
         system->configuration().particle(1).site(0).position()) < 1e-8,
    "The first site of each particle should be on top of each other.");
}

void Rotator::set_last_three_sites(const int ior, System * system) {
  ASSERT(system->configuration().domain().dimension() == 3, "hard coded for 3d");
  const Particle& part = system->configuration().particle(1);
  DEBUG("num sites:" << part.num_sites());
  DEBUG("ior: " << ior << " " << last_three_sites_.size());
  if (part.num_sites() >= 3) {
    const int num_sites = part.num_sites();
    for (int site = 0; site < 3; ++site) {
      const Position& site_pos = part.site(num_sites - 1 - site).position();
      for (int dim = 0; dim < 3; ++dim) {
        TRACE("ior " << ior << " site " << site << " dim " << dim);
        TRACE(last_three_sites_.size());
        TRACE(last_three_sites_[0].size());
        TRACE(last_three_sites_[0][0].dimension());
        last_three_sites_[ior][site].set_coord(dim, site_pos.coord(dim));
      }
    }
  } else if (part.num_sites() == 1) {
    for (int dim = 0; dim < 3; ++dim) {
      last_three_sites_[ior][0].set_coord(dim, part.site(0).position().coord(dim));
    }
  } else {
    FATAL("Particle has " << part.num_sites() << " sites but needs 1 or >= 3");
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
  if (tmp_3vec_[0] < 1e-8) {
    com1_.set_to_origin(3);
  } else {
    com1_.set_from_spherical(tmp_3vec_);
  }
  DEBUG("com1 " << com1_.str());
  select_->select_particle(1, system->configuration());
  rotate_->set_revert_possible(true, select_.get());
  if (eulers_.size() > 0) {
    DEBUG("euler " << eulers_[ior].str());
    eulers_[ior].compute_rotation_matrix(&rot_mat_);
    rotate_->move(*origin_, rot_mat_, system, select_.get());
  }
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

double Rotator::energy(const int ior, const double displacement, System * system) {
  update_xyz(ior, displacement, system);
  double en;
  ASSERT(select_->mobile().num_particles() == 1, "err");
  en = system->perturbed_energy(select_->mobile());
  revert(system);
  DEBUG("position of first site of mobile after revert " << system->configuration().particle(1).site(0).position().str());
  return en;
}

class ContactObjective : public Formula {
 public:
  ContactObjective(Rotator * rotator, System * system, const double ior, const double hard_u_limit, const double hard_u) {
    rotator_ = rotator;
    system_ = system;
    ior_ = ior;
    hard_u_limit_ = hard_u_limit;
    hard_u_ = hard_u;
  }

  double evaluate(const double distance) const override;
 private:
  Rotator * rotator_;
  System * system_;
  int ior_;
  double hard_u_limit_, hard_u_;
};

double Rotator::contact_distance(const int ior, System * system) {
  if (unique_[ior] != -1) {
    contact_[ior] = contact_[unique_[ior]];
    return contact_[unique_[ior]];
  }
  GoldenSearch minimize({{"tolerance", str(contact_tolerance_)},
    {"lower", "0"},
    {"upper", str(system->configuration().domain().max_side_length()/2)}});
  ContactObjective objective(this, system, ior, hard_limit_u_, hard_u_);
  const double dist = minimize.minimum(&objective) + 2.*contact_tolerance_;
  contact_[ior] = dist;
  if (!contact_xyz_file_name_.empty()) {
    update_xyz(ior, dist, system);
    contact_f_.write(contact_xyz_file_name_, system->configuration());
    revert(system);
  }
  //INFO("ior " << ior << " sph " << stheta_[ior] << " " << sphi_[ior] << " euler " << eulers_[ior].str());
  return dist;
}

double ContactObjective::evaluate(const double distance) const {
  const double en = rotator_->energy(ior_, distance, system_);
  TRACE("dist " << distance << " en " << en << " hu " << hard_u_limit_);
  if (hard_u_limit_ == hard_u_) {
    TRACE("hard contact");
    if (en > 1) {
      return 1e15/(distance + 0.1);
    } else {
      return distance + 0.1;
    }
  } else {
    TRACE("soft contact");
    if (en < hard_u_) {
      if (en > 0.) {
        const double obj = std::pow(en - hard_u_limit_, 2);
        TRACE("obj: " << obj);
        return obj;
      } else {
        const double obj = 1e101*distance;
        TRACE("obj: " << obj);
        return obj;
      }
    } else {
      return 1e101/(distance + 0.1);
    }
  }
}

void Rotator::gen_unique_orientations(const int num_orientations_per_pi, System * system, std::vector<std::vector<double> > bounds) {
  const Configuration& config = system->configuration();
  gen_orientations(num_orientations_per_pi, config, bounds);
  const double displacement = 40.;
  DEBUG("num orientations: " << num_orientations());
  ASSERT(num_proc_ == 1, "unique orientation search is not parallelized in the same way as the rest.");
  DEBUG("Set last three");
  INFO("num ori: " << num_orientations());
  //#pragma omp parallel for
  for (int ior = 0; ior < num_orientations(); ++ior) {
    update_xyz(ior, displacement, system);
    set_last_three_sites(ior, system);
    revert(system);
  }
  if (eulers_.size() > 0) {
    check_last_three_sites(0, system);
  }
  DEBUG("Determining unique orientations.");
//      //#pragma omp parallel for
  for (int ior = 0; ior < num_orientations(); ++ior) {
    unique_[ior] = -2;
  }
  std::vector<int> iors;
  #pragma omp parallel shared(iors)
  {
    auto thread = MakeThreadOMP();
    const int num_threads = thread->num();
    const int proc = thread->thread();
    std::unique_ptr<ProgressReport> report;
    if (proc == 0) {
      iors.resize(num_threads);
      report = std::make_unique<ProgressReport>(argtype({
        {"num", str(int(num_orientations()/num_threads))},
        {"task", "determine unique orientations"}}));
    }
    #pragma omp barrier
    for (int ior = proc; ior < num_orientations(); ior += num_threads) {
      determine_if_unique(ior, iors, num_threads, system);
      iors[proc] = ior;
      if (proc == 0) {
        report->check();
      }
    }
  }
  //#pragma omp parallel for schedule(static,1)
  //for (int ior = 0; ior < num_orientations(); ++ior) {
  //  determine_if_unique(ior, system);
  //}
}

}  // namespace feasst
