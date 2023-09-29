
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
#include "aniso/include/tabulate_two_rigid_body_3D.h"

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

TabulateTwoRigidBody3D::TabulateTwoRigidBody3D(argtype * args) {
  class_name_ = "TabulateTwoRigidBody3D";
  rotator_ = Rotator(args);
  num_orientations_per_pi_ = integer("num_orientations_per_pi", args, -1);
  num_z_ = integer("num_z", args, -1);
  ASSERT(num_z_ >= 2 || num_z_ == -1,
    "num_z:" << num_z_ << " must be >= 2 or ==-1");
  gamma_ = dble("gamma", args, -4);
  smoothing_distance_ = dble("smoothing_distance", args, 2);
  max_energy_ = dble("max_energy", args, 1e30);
  max_energy_set_ = dble("max_energy_set", args, 5);
  output_orientation_file_ = str("output_orientation_file", args, "");
  input_orientation_file_ = str("input_orientation_file", args, "");
  output_table_file_ = str("output_table_file", args, "");
  input_table_file_ = str("input_table_file", args, "");
  xyz_file_ = str("xyz_file", args, "");
  contact_xyz_file_ = str("contact_xyz_file", args, "");
}
TabulateTwoRigidBody3D::TabulateTwoRigidBody3D(argtype args) : TabulateTwoRigidBody3D(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapTabulateTwoRigidBody3D {
 public:
  MapTabulateTwoRigidBody3D() {
    auto obj = MakeTabulateTwoRigidBody3D({{"num_orientations_per_pi", "1"}});
    obj->deserialize_map()["TabulateTwoRigidBody3D"] = obj;
  }
};

static MapTabulateTwoRigidBody3D mapper_TabulateTwoRigidBody3D = MapTabulateTwoRigidBody3D();

TabulateTwoRigidBody3D::TabulateTwoRigidBody3D(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 8054 && version <= 8054, "mismatch version: " << version);
  feasst_deserialize(&num_orientations_per_pi_, istr);
  feasst_deserialize(&num_z_, istr);
  feasst_deserialize(&gamma_, istr);
  feasst_deserialize(&smoothing_distance_, istr);
  feasst_deserialize(&max_energy_, istr);
  feasst_deserialize(&max_energy_set_, istr);
  feasst_deserialize(&output_orientation_file_, istr);
  feasst_deserialize(&input_orientation_file_, istr);
  feasst_deserialize(&output_table_file_, istr);
  feasst_deserialize(&input_table_file_, istr);
  feasst_deserialize(&xyz_file_, istr);
  feasst_deserialize(&contact_xyz_file_, istr);
//  feasst_deserialize(&num_proc_, istr);
//  feasst_deserialize(&proc_, istr);
}

void TabulateTwoRigidBody3D::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(8054, ostr);
  feasst_serialize(num_orientations_per_pi_, ostr);
  feasst_serialize(num_z_, ostr);
  feasst_serialize(gamma_, ostr);
  feasst_serialize(smoothing_distance_, ostr);
  feasst_serialize(max_energy_, ostr);
  feasst_serialize(max_energy_set_, ostr);
  feasst_serialize(output_orientation_file_, ostr);
  feasst_serialize(input_orientation_file_, ostr);
  feasst_serialize(output_table_file_, ostr);
  feasst_serialize(input_table_file_, ostr);
  feasst_serialize(xyz_file_, ostr);
  feasst_serialize(contact_xyz_file_, ostr);
//  feasst_serialize(num_proc_, ostr);
//  feasst_serialize(proc_, ostr);
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

double TabulateTwoRigidBody3D::max_cubic_side_length(const int particle_type,
    const Configuration& config) const {
  double length = 0.;
  const Particle& part = config.particle_type(particle_type);
  for (const Site& site : part.sites()) {
    const double dist_sq = site.position().squared_distance();
    if (dist_sq > length*length) {
      length = std::sqrt(dist_sq);
    }
  }
  length *= 2; // double for +/- (e.g., radius -> diameter)
  const double max_sigma = config.model_params().select("sigma").mixed_max();
  const double max_cutoff = config.model_params().select("cutoff").mixed_max();
  //DEBUG("max_sigma " << max_sigma);
  //DEBUG("max_cutoff " << max_cutoff);
  length += 10 + max_sigma + max_cutoff;
  //length += 10 + max_sigma + max_sigma; // assuming hard sphere only
  length *= 2; // double for two particles.
  return length;
}

void TabulateTwoRigidBody3D::adjust_domain(System * system) {
  const double length0 = max_cubic_side_length(0, system->configuration());
  const double length1 = max_cubic_side_length(1, system->configuration());
  const double length = 0.5*(length0 + length1);
  const double delta_volume = std::pow(length, 3) - system->configuration().domain().volume();
  system->change_volume(delta_volume, {{"scale_particles", "false"}});
}

void TabulateTwoRigidBody3D::run(MonteCarlo * mc) {
  ASSERT(mc, "mc required");
  ASSERT(mc->configuration().dimension() == 3, "assumes 3D");
  System * system = mc->get_system();
  adjust_domain(system);
  rotator_.init(system, xyz_file_, contact_xyz_file_);
//  std::shared_ptr<ProgressReport> report;
  if (!input_table_file_.empty()) {
    read_contact_table(mc);
//    report = MakeProgressReport({{"num", str(rotator_.num_orientations())},
//                                 {"percent_per_write", "0.1"}});
  } else {
    DEBUG("Determining orientations.");
    if (input_orientation_file_.empty()) {
      rotator_.gen_orientations(num_orientations_per_pi_, mc->configuration());
//      report = MakeProgressReport({{"num", str(rotator_.num_orientations())},
//                                   {"percent_per_write", "0.1"}});
      const double displacement = 40.;

      DEBUG("num orientations: " << rotator_.num_orientations());
      ASSERT(rotator_.num_proc_ == 1, "unique orientation search is not parallelized in the same way as the rest.");
      DEBUG("Set last three");
      //#pragma omp parallel for
      for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
        rotator_.update_xyz(ior, displacement, system);
        rotator_.set_last_three_sites(ior, system);
        //DEBUG("ior " << ior << " unique " << rotator_.unique_[ior]); // << " last three " << last_three[0].str() << " " << last_three[1].str() << " " << last_three[2].str());
        rotator_.revert(system);
//        report->check();
      }
      rotator_.check_last_three_sites(0, system);
      DEBUG("Determining unique orientations.");
//      //#pragma omp parallel for
      for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
        rotator_.unique_[ior] = -2;
      }
      std::vector<int> iors;
      #pragma omp parallel shared(iors)
      {
        auto thread = MakeThreadOMP();
        const int num_threads = thread->num();
        const int proc = thread->thread();
        if (proc == 0) {
          iors.resize(num_threads);
        }
        #pragma omp barrier
        for (int ior = proc; ior < rotator_.num_orientations(); ior += num_threads) {
          rotator_.determine_if_unique(ior, iors, num_threads, system);
          iors[proc] = ior;
        }
      }
      //#pragma omp parallel for schedule(static,1)
      //for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
      //  rotator_.determine_if_unique(ior, system);
      //}
      if (!output_orientation_file_.empty()) {
        ouput_orientations_();
        DEBUG("num orientations: " << rotator_.num_orientations());
        DEBUG("num unique: " << rotator_.num_unique());
        DEBUG("fraction unique: " << rotator_.fraction_unique());
        return;
      }
//      report->reset();
    } else {
      read_input_orientations_(mc->configuration());
    }

    DEBUG("Obtaining contact distances.");
    for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
      rotator_.contact_distance(ior, system);
      //const double dist = rotator_.contact_distance(ior, system);
      //DEBUG("ior " << ior << " dist " << MAX_PRECISION << dist);
//      report->check();
    }
  }

  if (num_z_ != -1) {
    DEBUG("Obtaining energies for each orientation.");
    resize(rotator_.num_orientations(), num_z_, &rotator_.energy_);
    //report->reset();
    const double dz = 1./static_cast<double>(num_z_ - 1);
    const double rc = mc->configuration().model_params().select("cutoff").mixed_max();
    //for (int ior = 0; ior < 1; ++ior) {
    for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
      const int unique_ior = rotator_.unique_[ior];
      if (unique_ior == -1) {
        const double rh = rotator_.contact_[ior];
        const double rhg = std::pow(rh, gamma_);
        const double rcg = std::pow(rh + rc - smoothing_distance_, gamma_);
        for (int iz = 0; iz < num_z_; ++iz) {
          const double z = iz*dz;
          const double dist = std::pow(z*(rcg - rhg) + rhg, 1./gamma_);
          double en = rotator_.energy(ior, dist, system);
          if (en > max_energy_) {
            ASSERT(iz > 0, "ior " << ior << " z " << z << " dist " << dist <<
              " en " << en << " . Incorrect contact distance?");
            en = max_energy_set_;
          }
          rotator_.energy_[ior][iz] = en;
        }
        //const double dist = rotator_.contact_distance(ior, system);
        //DEBUG("ior " << ior);
//      } else {
//        for (int iz = 0; iz < num_z_; ++iz) {
//          rotator_.energy_[ior][iz] = rotator_.energy_[unique_ior][iz];
//        }
      }
//        report->check();
    }
  }
  DEBUG("Outputing table");
  write_table(mc);
}

void TabulateTwoRigidBody3D::write_table(MonteCarlo * mc) const {
  ASSERT(!output_table_file_.empty(), "no output_table_file");
  std::ofstream file(output_table_file_);
  ASSERT(file.good(), "err");
  file << "site_types 1 0" << std::endl;
  std::streamsize ss = std::cout.precision();
  file << MAX_PRECISION << "num_orientations_per_pi " << num_orientations_per_pi_ << std::endl;
  file << "gamma " << gamma_ << std::endl;
  file << "delta ";
  if (num_z_ == -1) {
    file << "0" << std::endl;
  } else {
    file << mc->configuration().model_params().select("cutoff").mixed_max() << std::endl;
  }
  file << "num_z " << num_z_ << std::endl;
  file << "smoothing_distance " << smoothing_distance_ << std::endl;
  file << std::setprecision(ss);
  for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
    const int unique_ior = rotator_.unique_[ior];
    if (unique_ior == -1) {
      if (num_z_ == -1) {
        file << std::setprecision(std::numeric_limits<float>::digits10) << rotator_.contact_[ior];
      } else {
        file << std::scientific << rotator_.contact_[ior];
        for (int iz = 0; iz < num_z_; ++iz) {
          file << " " << rotator_.energy_[ior][iz];
        }
      }
    } else {
      file << "-1 " << unique_ior;
    }
    file << std::endl;
  }
}

void TabulateTwoRigidBody3D::read_contact_table(MonteCarlo * mc) {
  std::ifstream file(input_table_file_);
  ASSERT(file.good(), "Error reading file: " << input_table_file_);
  std::string line;
  std::getline(file, line);
  DEBUG("line " << line);
  file >> line >> num_orientations_per_pi_;
  DEBUG("num_orientations_per_pi_ " << num_orientations_per_pi_);
  ASSERT(num_orientations_per_pi_ >= 1,
    "num_orientations_per_pi:" << num_orientations_per_pi_);
  rotator_.gen_orientations(num_orientations_per_pi_, mc->configuration());
  for (int i = 0; i < 5; ++i) {
    std::getline(file, line);
    DEBUG("line " << i << ": " << line);
  }
  int ior = 0;
  for (int iorall = 0; iorall < rotator_.num_orientations_all_proc(); ++iorall) {
    double first_value;
    int unique_ior;
    file >> first_value;
    if (first_value < 0) {
      file >> unique_ior;
    }
    if (rotator_.ior_in_proc(iorall)) {
      if (first_value < 0) {
        //rotator_.contact_[ior] = rotator_.contact_[unique_ior];
        rotator_.unique_[ior] = unique_ior;
      } else {
        rotator_.contact_[ior] = first_value;
        rotator_.unique_[ior] = -1;
      }
    //INFO("rotator_.unique_[ior] " << rotator_.unique_[ior]);
      ++ior;
    }
  }
}

void TabulateTwoRigidBody3D::ouput_orientations_() {
  std::ofstream file(output_orientation_file_);
  ASSERT(file.good(), "err");
  file << "num_orientations_per_pi " << num_orientations_per_pi_ << std::endl;
  for (int ior = 0; ior < rotator_.num_orientations(); ++ior) {
    file << rotator_.unique_[ior] << " ";
  }
}

void TabulateTwoRigidBody3D::read_input_orientations_(const Configuration& config) {
  std::ifstream file(input_orientation_file_);
  ASSERT(file.good(), "Cannot read " << input_orientation_file_);
  int num_o_p_pi;
  std::string tmp_str;
  file >> tmp_str >> num_o_p_pi;
  if (num_orientations_per_pi_ == -1) {
    num_orientations_per_pi_ = num_o_p_pi;
  } else {
    ASSERT(num_orientations_per_pi_ == num_o_p_pi,
      "The given num_orientations_per_pi: " << num_orientations_per_pi_ <<
      " is not equal to the one in the file: " << num_o_p_pi);
  }
  rotator_.gen_orientations(num_orientations_per_pi_, config);
  //INFO(rotator_.num_orientations());
  ASSERT(rotator_.num_orientations() > 0, "err");
  int ior = 0;
  for (int iorall = 0; iorall < rotator_.num_orientations_all_proc(); ++iorall) {
    int unique;
    file >> unique;
    if (rotator_.ior_in_proc(iorall)) {
      rotator_.unique_[ior] = unique;
      ASSERT(file.peek() != EOF, "orientation file has less orientations than "
        << "expected. Is num_orientations_per_pi the same?");
      //INFO(rotator_.unique_[ior]);
      ++ior;
    }
  }
  if (rotator_.num_proc() == 1) {
    tmp_str.clear();
    file >> tmp_str;
    ASSERT(tmp_str.empty() && file.peek() == EOF,
      "orientation file has more orientations than expected. " <<
      "Is num_orientations_per_pi the same? " <<
      "Last line was read as: " << tmp_str);
  }
}

void TabulateTwoRigidBody3D::combine_table(argtype args) const {
  const int num_header = 6;
  const int num = rotator_.num_proc();
  const std::string prefix = str("prefix", args);
  const std::string suffix = str("suffix", args);
  std::ofstream output(prefix + suffix);
  std::vector<std::ifstream> input(num);
  std::string line;
  for (int proc = 0; proc < num; ++proc) {
    input[proc] = std::ifstream(prefix + str(proc) + suffix);
    for (int i = 0; i < num_header; ++i) {
      std::getline(input[proc], line);
      if (proc == 0) {
        output << line << std::endl;
      }
    }
  }
  bool done = false;
  int proc = 0;
  int i = 0;
  while (!done) {
    std::getline(input[proc], line);
    output << line << std::endl;
    ++proc;
    if (proc == num) {
      proc = 0;
    }
    ASSERT(i < 1e30, "infinite loop?");
    ++i;
    if (i == rotator_.num_orientations_all_proc()) {
      done = true;
    }
  }
}

}  // namespace feasst
