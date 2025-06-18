
#include <iostream>
#include <chrono> // sleep
#include <thread> // sleep
#include "utils/include/serialize.h"
#include "utils/include/utils.h"
#include "utils/include/arguments_extra.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/formula.h"
#include "math/include/golden_search.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/rotator.h"
#include "aniso/include/tabulate_two_rigid_body_3D.h"

namespace feasst {

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
  contact_xyz_index_ = integer("contact_xyz_index", args, -1);
}
TabulateTwoRigidBody3D::TabulateTwoRigidBody3D(argtype args) : TabulateTwoRigidBody3D(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(TabulateTwoRigidBody3D, argtype({{"num_orientations_per_pi", "1"}}));

TabulateTwoRigidBody3D::TabulateTwoRigidBody3D(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 8054 && version <= 8055, "mismatch version: " << version);
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
  if (version >= 8055) {
    feasst_deserialize(&contact_xyz_index_, istr);
  }
//  feasst_deserialize(&num_proc_, istr);
//  feasst_deserialize(&proc_, istr);
}

void TabulateTwoRigidBody3D::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(8055, ostr);
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
  feasst_serialize(contact_xyz_index_, ostr);
//  feasst_serialize(num_proc_, ostr);
//  feasst_serialize(proc_, ostr);
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
  DEBUG("length " << length);
  const double volume = system->configuration().domain().volume();
  DEBUG("volume " << volume);
  const double delta_volume = std::pow(length, 3) - volume;
  DEBUG("delta_volume " << delta_volume);
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
    int ior_first = 0;
    int ior_less_than = rotator_.num_orientations();
    if (contact_xyz_index_ != -1) {
      ior_first = contact_xyz_index_;
      ior_less_than = ior_first + 1;
    }
    for (int ior = ior_first; ior < ior_less_than; ++ior) {
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
    input[proc].open(prefix + str(proc) + suffix);
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
