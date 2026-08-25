#include <fstream>
#include "threads/include/thread_omp.h"
#include "utils/include/file.h"  // num_lines
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "math/include/accumulator.h"
#include "math/include/constants.h"
#include "math/include/recursive_table.h"
#include "math/include/utils_math.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/rotator.h"
#include "aniso/include/recursive_table.h"
#include "aniso/include/model_recursive_table.h"
#include "aniso/include/build_recursive_table.h"

namespace feasst {

BuildRecursiveTable::BuildRecursiveTable(argtype * args) {
  class_name_ = "BuildRecursiveTable";
  mayer_training_file_ = str("mayer_training_file", args, "");
  cutoff_ = dble("cutoff", args, -1);
  output_file_ = str("output_file", args, "");
  verbose_file_ = str("verbose_file", args, "");
  extra_verbose_ = boolean("extra_verbose", args, false);
  hard_limit_u_ = dble("hard_limit_u", args, 100);
  num_z_ = dble("num_z", args, 5);
  num_orientations_per_pi_ = dble("num_orientations_per_pi", args, 5);
  beta_ = dble("beta", args, 1.);
  min_criteria_ = dble("min_criteria", args, 0.03);
  min_criteria_energy_ = dble("min_criteria_energy", args, 0.);
  if (min_criteria_energy_ <= 0.) {
    min_criteria_energy_ = min_criteria_;
  }
  contact_only_ = boolean("contact_only", args, false);
  energy_only_ = boolean("energy_only", args, false);
  input_recursive_table_ = str("input_recursive_table", args, "");
  assume_all_unique_ = str("assume_all_unique", args, "false");
  processor_ = integer("processor", args, 0);
  num_processors_ = integer("num_processors", args, 1);
  base_table_file_ = str("base_table_file", args, "");
  trained_file_ = str("trained_file", args, "");
  stage_ = str("stage", args, "all");
  if (stage_ == "base") {
    min_criteria_ = NEAR_INFINITY_FLOAT;
    min_criteria_energy_ = NEAR_INFINITY_FLOAT;
    assume_all_unique_ = "true";
  } else if (stage_.substr(0,5) == "build") {
    INFO("setting assume_all_unique true");
    assume_all_unique_ = "true";
  }
}
BuildRecursiveTable::BuildRecursiveTable(argtype args) : BuildRecursiveTable(&args) {
  feasst_check_all_used(args);
}
BuildRecursiveTable::~BuildRecursiveTable() {}

FEASST_MAPPER(BuildRecursiveTable,);

BuildRecursiveTable::BuildRecursiveTable(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3257 && version <= 3257, "mismatch version: " << version);
  feasst_deserialize(&mayer_training_file_, istr);
  feasst_deserialize(&cutoff_, istr);
  feasst_deserialize(&output_file_, istr);
  feasst_deserialize(&verbose_file_, istr);
  feasst_deserialize(&hard_limit_u_, istr);
  feasst_deserialize(&num_z_, istr);
  feasst_deserialize(&num_orientations_per_pi_, istr);
  feasst_deserialize(&beta_, istr);
  feasst_deserialize(&min_criteria_, istr);
  feasst_deserialize(&min_criteria_energy_, istr);
  feasst_deserialize(&contact_only_, istr);
  feasst_deserialize(&energy_only_, istr);
  feasst_deserialize(&extra_verbose_, istr);
  feasst_deserialize(&assume_all_unique_, istr);
  feasst_deserialize(&processor_, istr);
  feasst_deserialize(&num_processors_, istr);
  feasst_deserialize(&base_table_file_, istr);
  feasst_deserialize(&trained_file_, istr);
  feasst_deserialize(&stage_, istr);
}

void BuildRecursiveTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_build_recursive_table_(ostr);
}

void BuildRecursiveTable::serialize_build_recursive_table_(std::ostream& ostr) const {
  serialize_action_(ostr);
  feasst_serialize_version(3257, ostr);
  feasst_serialize(mayer_training_file_, ostr);
  feasst_serialize(cutoff_, ostr);
  feasst_serialize(output_file_, ostr);
  feasst_serialize(verbose_file_, ostr);
  feasst_serialize(hard_limit_u_, ostr);
  feasst_serialize(num_z_, ostr);
  feasst_serialize(num_orientations_per_pi_, ostr);
  feasst_serialize(beta_, ostr);
  feasst_serialize(min_criteria_, ostr);
  feasst_serialize(min_criteria_energy_, ostr);
  feasst_serialize(contact_only_, ostr);
  feasst_serialize(energy_only_, ostr);
  feasst_serialize(extra_verbose_, ostr);
  feasst_serialize(assume_all_unique_, ostr);
  feasst_serialize(processor_, ostr);
  feasst_serialize(num_processors_, ostr);
  feasst_serialize(base_table_file_, ostr);
  feasst_serialize(trained_file_, ostr);
  feasst_serialize(stage_, ostr);
}

void BuildRecursiveTable::build_table_(const std::vector<std::vector<double> >& bounds, RecursiveTable1D * tab, MonteCarlo * mc) const {
  INFO("building table. inner:" << bounds[0][0] << " outer:" << bounds[0][1]);
  ASSERT(static_cast<int>(bounds.size()) == 1, "mismatch");
  const double lower = bounds[0][0];
  const double upper = bounds[0][1];
  // store the positions of the two single-site particles
  std::vector<std::vector<double> > coords(2), old_coords(2);
  const Configuration& config = mc->system().configuration();
  coords[0] = {0, 0, 0};
  old_coords[0] = config.particle(0).site(0).position().coord();
  coords[1] = {0, 0, 0};
  old_coords[1] = config.particle(1).site(0).position().coord();
  ASSERT(config.particle(0).site(0).position().squared_distance() < 1e-8, "Err");
  std::vector<std::vector<double> > coords_old = coords;
  System * sys = mc->get_system();
  for (int i = 0; i < tab->num(); ++i) {
    if (i % num_processors_ == processor_) {
      const double z = static_cast<double>(i)/(tab->num() - 1);
      const double r = lower + z*(upper - lower);
      DEBUG("r:" << r);
      coords[1][0] = r;
      sys->get_configuration()->update_positions(coords);
      const double en = sys->energy();
      tab->set_data(i, en);
    }
  }
  sys->get_configuration()->update_positions(old_coords);
  sys->energy();
  DEBUG("table built");
}

std::string BuildRecursiveTable::verbose_name_(const int iteration) const {
  if (verbose_file_.empty()) {
    return std::string("");
  }
  return verbose_file_ + "_" + str(iteration);
}

int energy_index(const std::vector<std::vector<double> >& data) {
  DEBUG(data[0].size());
  if (data[0].size() >= 7) {
    return 6;
  }
  return static_cast<int>(data[0].size()) - 1;
}

RecursiveTable1D BuildRecursiveTable::build_1dcontact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system, std::vector<System> * systems, const bool cutoff) {
  Rotator rotator = gen_rotator_(bounds, system);
  INFO("num orientations:" << rotator.num_orientations());
  int num_threads;
  std::vector<Rotator> rotators;
  gen_rotators_(bounds, systems, &num_threads, &rotators);
  RecursiveTable1D contact(argtype({{"num", str(rotator.sizes_[0])}}));
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf);
    file << "r,theta,x,y" << std::endl;
  }
  #pragma omp parallel shared(contact)
  {
    Position pos(argtype({{"dimension", "2"}}));
    const int ithread = ThreadOMP().thread();
    if (ithread < num_threads) {
      Rotator * irotator = get_rotator_(ithread, &rotator, &rotators);
      System * isystem = get_system_(ithread, system, systems);
      for (int ior = ithread; ior < irotator->num_orientations(); ior += num_threads) {
        if (ior % num_processors_*num_threads == processor_) {
          double dist;
          if (cutoff) {
            dist = irotator->cutoff_distance(ior, isystem);
          } else{
            dist = irotator->contact_distance(ior, isystem);
          }
          const std::vector<int>& idx = irotator->indices_[ior];
          contact.set_data(idx[0], dist);
          if (!verbf.empty()) {
            #pragma omp critical
            {
              pos.set_from_spherical({dist, irotator->stheta_[ior]});
              file << dist << "," << irotator->stheta_[ior] << "," << pos.str()
                   << std::endl;
            }
          }
        }
      }
    }
  }
  return contact;
}

RecursiveTable2D BuildRecursiveTable::build_2dcontact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system, std::vector<System> * systems, const bool cutoff) {
  Rotator rotator = gen_rotator_(bounds, system);
  INFO("num orientations:" << rotator.num_orientations());
  int num_threads;
  std::vector<Rotator> rotators;
  gen_rotators_(bounds, systems, &num_threads, &rotators);
  RecursiveTable2D contact(argtype({{"num0", str(rotator.sizes_[0])},
                                    {"num1", str(rotator.sizes_[1])}}));
  std::ofstream file;
  const int dimen = system->configuration().dimension();
  if (!verbf.empty()) {
    file.open(verbf);
    if (dimen == 3) {
      file << "r,theta,phi,x,y,z" << std::endl;
    } else {
      file << "r,theta,phi,x,y" << std::endl;
    }
  }
  #pragma omp parallel shared(contact)
  {
    Position pos(argtype({{"dimension", str(dimen)}}));
    const int ithread = ThreadOMP().thread();
    if (ithread < num_threads) {
      Rotator * irotator = get_rotator_(ithread, &rotator, &rotators);
      System * isystem = get_system_(ithread, system, systems);
      for (int ior = ithread; ior < irotator->num_orientations(); ior += num_threads) {
        if (ior % num_processors_*num_threads == processor_) {
          double dist;
          if (cutoff) {
            dist = irotator->cutoff_distance(ior, isystem);
          } else{
            dist = irotator->contact_distance(ior, isystem);
          }
          DEBUG("dist " << dist);
          const std::vector<int>& idx = irotator->indices_[ior];
          DEBUG("idx " << feasst_str(idx));
          contact.set_data(idx[0], idx[1], dist);
          if (!verbf.empty()) {
            #pragma omp critical
            {
              if (dimen == 3) {
                pos.set_from_spherical({dist, irotator->stheta_[ior], irotator->sphi_[ior]});
                file << dist << "," << irotator->stheta_[ior] << "," << irotator->sphi_[ior]
                  << "," << pos.str() << std::endl;
              } else {
                pos.set_from_spherical({dist, irotator->stheta_[ior]});
                file << dist << "," << irotator->stheta_[ior] << "," << irotator->eulers_[ior].phi()
                  << "," << pos.str() << std::endl;
              }
            }
          }
        }
      }
    }
  }
  return contact;
}

// HWH copied from below
RecursiveTable6D BuildRecursiveTable::build_energy_(const std::vector<std::vector<double> >& bounds, const std::vector<double>& zbnds, const std::string& verbf, const RecursiveTable5D& contact, const RecursiveTable5D& cutoff, System * system, std::vector<System> * systems) {
  Rotator rotator = gen_rotator_(bounds, system);
  INFO("num orientations:" << rotator.num_orientations());
  int num_threads;
  std::vector<Rotator> rotators;
  gen_rotators_(bounds, systems, &num_threads, &rotators);
  RecursiveTable6D en(argtype({{"num0", str(rotator.sizes_[0])},
                               {"num1", str(rotator.sizes_[1])},
                               {"num2", str(rotator.sizes_[2])},
                               {"num3", str(rotator.sizes_[3])},
                               {"num4", str(rotator.sizes_[4])},
                               {"num5", str(num_z_)}}));
  std::ofstream file;
  const int dimen = system->configuration().dimension();
  if (!verbf.empty()) {
    file.open(verbf);
    file << "r,theta,phi,ephi,etheta,epsi,x,y,z,en" << std::endl;
  }
  #pragma omp parallel shared(en)
  {
    Position pos(argtype({{"dimension", str(dimen)}}));
    double s1, s2, s3, s4, s5;
    const double dzpn = (zbnds[1] - zbnds[0])/static_cast<double>(num_z_ - 1);
    bool full_range = false;
    int subsize = num_z_;
    int idis_start = 0;
    if (std::abs(zbnds[1] - zbnds[0] - 1) < 1e-8) {
      full_range = true;
      --subsize;
      idis_start = 1;
    }
    INFO("full range " << full_range << " subsize " << subsize << " dzpn " << dzpn);
    const int ithread = ThreadOMP().thread();
    if (ithread < num_threads) {
      Rotator * irotator = get_rotator_(ithread, &rotator, &rotators);
      System * isystem = get_system_(ithread, system, systems);
      for (int ior = ithread; ior < irotator->num_orientations(); ior += num_threads) {
        const std::vector<int>& idx = irotator->indices_[ior];
        const Euler& eul = irotator->eulers_[ior];
        scaled_relative_orientation(irotator->stheta_[ior], irotator->sphi_[ior], eul.phi(), eul.theta(), eul.psi(), &s1, &s2, &s3, &s4, &s5);
        //INFO("s1 " << s1 << " s2 " << s2);
        const double rh = contact.linear_interpolation(s1, s2, s3, s4, s5);
        const double rc = cutoff.linear_interpolation(s1, s2, s3, s4, s5);
        //INFO("rh " << rh);
        //INFO("rc " << rc);
        for (int idis = idis_start; idis < subsize; ++idis) {
          DEBUG("idis " << idis << " ior " << ior << " num_z " << num_z_ << " idis + ior*numz " << (idis + ior*num_z_) << " np*nt " << (num_processors_*num_threads) << " con " << (idis + ior*num_z_) % (num_processors_*num_threads) << " proc " << processor_);
          if ( (num_threads != 1) || (idis + ior*num_z_) % (num_processors_*num_threads) == processor_) {
            const double z = zbnds[0] + dzpn*static_cast<double>(idis);
            const double dist = rh + z*(rc - rh);
            double ener = irotator->energy(ior, dist, isystem);
            if (ener > NEAR_INFINITY_FLOAT) {
              ener = NEAR_INFINITY_FLOAT;
            }
            en.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], idis, ener);
            if (!verbf.empty()) {
              #pragma omp critical
              {
                pos.set_from_spherical({dist, irotator->stheta_[ior], irotator->sphi_[ior]});
                file << dist << "," << irotator->stheta_[ior] << "," << irotator->sphi_[ior]
                  << "," << eul.phi() << "," << eul.theta() << "," << eul.psi() << ","
                  << pos.str() << "," << ener << std::endl;
              }
            }
          }
        }
        if (full_range) {
          /* Set energy to zero at contact and cutoff (z=0 and 1)
             Although it may seem more efficient to not even store these known
             values, these points are non trivial and required in the nested
             tables. If the majority of data points are nested then the
             potential efficiency gain is minor, and more book keeping would be
             required to interpolate differently from base vs non-nested data
             points. The energy is not computed so tabling is not slower this
             way, only the data table size will be a bit larger.
           */
          en.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], 0, hard_limit_u_);
          en.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], num_z_ - 1, 0.);
        }
      }
    }
  }
  return en;
}

RecursiveTable3D BuildRecursiveTable::build_3denergy_(const std::vector<std::vector<double> >& bounds, const std::vector<double>& zbnds, const std::string& verbf, const RecursiveTable2D& contact, const RecursiveTable2D& cutoff, System * system, std::vector<System> * systems) {
  Rotator rotator = gen_rotator_(bounds, system);
  INFO("num orientations:" << rotator.num_orientations());
  int num_threads;
  std::vector<Rotator> rotators;
  gen_rotators_(bounds, systems, &num_threads, &rotators);
  RecursiveTable3D en(argtype({{"num0", str(rotator.sizes_[0])},
                               {"num1", str(rotator.sizes_[1])},
                               {"num2", str(num_z_)}}));
  //ASSERT(contact.percent_nested() > 0, "testing?");
  std::ofstream file;
  const int dimen = system->configuration().dimension();
  if (!verbf.empty()) {
    file.open(verbf);
    if (dimen == 3) {
      file << "r,theta,phi,x,y,z,en" << std::endl;
    } else {
      file << "r,theta,phi,x,y,en" << std::endl;
    }
  }
  #pragma omp parallel shared(en)
  {
    Position pos(argtype({{"dimension", str(dimen)}}));
    double s1, s2;
    const double dzpn = (zbnds[1] - zbnds[0])/static_cast<double>(num_z_ - 1);
    bool full_range = false;
    int subsize = num_z_;
    int idis_start = 0;
    if (std::abs(zbnds[1] - zbnds[0] - 1) < 1e-8) {
      full_range = true;
      --subsize;
      idis_start = 1;
    }
    DEBUG("full range " << full_range << " subsize " << subsize << " dzpn " << dzpn);
    const int ithread = ThreadOMP().thread();
    if (ithread < num_threads) {
      Rotator * irotator = get_rotator_(ithread, &rotator, &rotators);
      System * isystem = get_system_(ithread, system, systems);
      for (int ior = ithread; ior < irotator->num_orientations(); ior += num_threads) {
        const std::vector<int>& idx = irotator->indices_[ior];
        if (dimen == 3) {
          scaled_relative_orientation(irotator->stheta_[ior], irotator->sphi_[ior], dimen, &s1, &s2);
        } else if (dimen == 2) {
          scaled_relative_orientation(irotator->stheta_[ior], irotator->eulers_[ior].phi(), dimen, &s1, &s2);
        } else {
          FATAL("unrecognized dimen");
        }
        const double rh = contact.linear_interpolation(s1, s2);
        const double rc = cutoff.linear_interpolation(s1, s2);
        //for (int idis = 0; idis < subsize; ++idis) {
        for (int idis = idis_start; idis < subsize; ++idis) {
          //if (ior % num_processors_*num_threads == processor_) {
          if ( (num_threads != 1) || (idis + ior*num_z_) % (num_processors_*num_threads) == processor_) {
            const double z = zbnds[0] + dzpn*static_cast<double>(idis);
            const double dist = rh + z*(rc - rh);
            double ener = irotator->energy(ior, dist, isystem);
            if (ener > NEAR_INFINITY_FLOAT) {
              ener = NEAR_INFINITY_FLOAT;
            }
            en.set_data(idx[0], idx[1], idis, ener);
            if (!verbf.empty()) {
              #pragma omp critical
              {
                file << dist << "," << irotator->stheta_[ior];
                if (dimen == 3) {
                  pos.set_from_spherical({dist, irotator->stheta_[ior], irotator->sphi_[ior]});
                  file << "," << irotator->sphi_[ior];
                } else if (dimen == 2) {
                  pos.set_from_spherical({dist, irotator->stheta_[ior]});
                  file << "," << irotator->eulers_[ior].phi();
                } else {
                  FATAL("unrecognized dimen");
                }
                file << "," << pos.str() << "," << ener << std::endl;
              }
            }
          }
          if (full_range) {
            // Set energy to zero at cutoff (z=1)
            en.set_data(idx[0], idx[1], 0, hard_limit_u_);
            en.set_data(idx[0], idx[1], num_z_ - 1, 0.);
          }
        }
      }
    }
  }
  return en;
}

// HWH copied from above
RecursiveTable2D BuildRecursiveTable::build_2denergy_(const std::vector<std::vector<double> >& bounds, const std::vector<double>& zbnds, const std::string& verbf, const RecursiveTable1D& contact, const RecursiveTable1D& cutoff, System * system, std::vector<System> * systems) {
  Rotator rotator = gen_rotator_(bounds, system);
  INFO("num orientations:" << rotator.num_orientations());
  int num_threads;
  std::vector<Rotator> rotators;
  gen_rotators_(bounds, systems, &num_threads, &rotators);
  RecursiveTable2D en(argtype({{"num0", str(rotator.sizes_[0])},
                               {"num1", str(num_z_)}}));
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf);
    file << "r,theta,x,y,en" << std::endl;
  }
  #pragma omp parallel shared(en)
  {
    Position pos(argtype({{"dimension", "2"}}));
    INFO("num " << rotator.num_orientations()*num_z_);
    double s1;
    const double dzpn = (zbnds[1] - zbnds[0])/static_cast<double>(num_z_ - 1);
    bool full_range = false;
    int subsize = num_z_;
    int idis_start = 0;
    if (std::abs(zbnds[1] - zbnds[0] - 1) < 1e-8) {
      full_range = true;
      --subsize;
      idis_start = 1;
    }
    INFO("full range " << full_range << " subsize " << subsize << " dzpn " << dzpn);
    const int ithread = ThreadOMP().thread();
    if (ithread < num_threads) {
      Rotator * irotator = get_rotator_(ithread, &rotator, &rotators);
      System * isystem = get_system_(ithread, system, systems);
      for (int ior = ithread; ior < irotator->num_orientations(); ior += num_threads) {
        const std::vector<int>& idx = irotator->indices_[ior];
        DEBUG("stheta " << irotator->stheta_[ior]);
        scaled_relative_orientation(irotator->stheta_[ior], &s1);
        DEBUG("s1 " << s1);
        const double rh = contact.linear_interpolation(s1);
        //INFO("rh " << rh);
        const double rc = cutoff.linear_interpolation(s1);
        //INFO("rc " << rc);
        for (int idis = idis_start; idis < subsize; ++idis) {
          if ( (num_threads != 1) || (idis + ior*num_z_) % (num_processors_*num_threads) == processor_) {
            const double z = zbnds[0] + dzpn*static_cast<double>(idis);
            const double dist = rh + z*(rc - rh);
            double ener = irotator->energy(ior, dist, isystem);
            if (ener > NEAR_INFINITY_FLOAT) {
              ener = NEAR_INFINITY_FLOAT;
            }
            en.set_data(idx[0], idis, ener);
            if (!verbf.empty()) {
              #pragma omp critical
              {
                pos.set_from_spherical({dist, irotator->stheta_[ior]});
                file << dist << "," << irotator->stheta_[ior]
                  << "," << pos.str() << "," << ener << std::endl;
              }
            }
          }
          if (full_range) {
            // Set energy to zero at cutoff (z=1)
            en.set_data(idx[0], 0, hard_limit_u_);
            en.set_data(idx[0], num_z_ - 1, 0.);
            if (!verbf.empty()) {
              #pragma omp critical
              {
                const int ior = subsize;
                scaled_relative_orientation(irotator->stheta_[ior], &s1);
                const double dist = cutoff.linear_interpolation(s1);
                pos.set_from_spherical({dist, irotator->stheta_[ior]});
                file << dist << "," << irotator->stheta_[subsize]
                  << "," << pos.str() << "," << 0. << std::endl;
              }
            }
          }
        }
      }
    }
  }
  return en;
}

Rotator BuildRecursiveTable::gen_rotator_(
    const std::vector<std::vector<double> >& bounds, System * system) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}, {"assume_all_unique", assume_all_unique_}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
  return rotator;
}

int BuildRecursiveTable::get_num_threads_() {
  int num_threads = 1;
  #pragma omp parallel
  {
    num_threads = ThreadOMP().num();
  }
  if (num_threads > 1 && num_processors_ > 1) {
    //WARN("Simultaneous parallelization with OMP and by processor is not implemented. Parallelizing by processor instead of OMP by setting num_threads to 1.");
    num_threads = 1;
  }
  return num_threads;
}

void BuildRecursiveTable::gen_rotators_(
    const std::vector<std::vector<double> >& bounds,
    std::vector<System> * systems,
    int * num_threads, std::vector<Rotator> * rotators) {
  *num_threads = get_num_threads_();
  rotators->resize(*num_threads - 1);
  ASSERT(systems->size() == rotators->size(), "mismatch sizes of system:"
      << systems->size() << " rotators:" << rotators->size());
  for (int thread = 0; thread < *num_threads - 1; ++thread) {
    (*rotators)[thread] = gen_rotator_(bounds, &(*systems)[thread]);
  }
}

Rotator * BuildRecursiveTable::get_rotator_(const int ithread, Rotator * rotator,
    std::vector<Rotator> * rotators) {
  Rotator * irotator;
  if (ithread == 0) {
    irotator = rotator;
  } else {
    irotator = &(*rotators)[ithread - 1];
  }
  ASSERT(irotator, "er");
  return irotator;
}

System * BuildRecursiveTable::get_system_(const int ithread, System * system,
    std::vector<System> * systems) {
  System * isystem;
  if (ithread == 0) {
    isystem = system;
  } else {
    isystem = &(*systems)[ithread - 1];
  }
  ASSERT(isystem, "er");
  return isystem;
}

RecursiveTable5D BuildRecursiveTable::build_contact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system, std::vector<System> * systems, const bool cutoff) {
  Rotator rotator = gen_rotator_(bounds, system);
  DEBUG("num orientations:" << rotator.num_orientations());
  int num_threads;
  std::vector<Rotator> rotators;
  gen_rotators_(bounds, systems, &num_threads, &rotators);
  RecursiveTable5D contact(argtype({{"num0", str(rotator.sizes_[0])},
                                    {"num1", str(rotator.sizes_[1])},
                                    {"num2", str(rotator.sizes_[2])},
                                    {"num3", str(rotator.sizes_[3])},
                                    {"num4", str(rotator.sizes_[4])}}));
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf);
    file << "r,theta,phi,ephi,etheta,epsi,x,y,z" << std::endl;
  }
  #pragma omp parallel shared(contact)
  {
    Position pos(argtype({{"dimension", "3"}}));
    const int ithread = ThreadOMP().thread();
    if (ithread < num_threads) {
      Rotator * irotator = get_rotator_(ithread, &rotator, &rotators);
      System * isystem = get_system_(ithread, system, systems);
      for (int ior = ithread; ior < irotator->num_orientations(); ior += num_threads) {
        if (ior % num_processors_*num_threads == processor_) {
          double dist;
          if (cutoff) {
            dist = irotator->cutoff_distance(ior, isystem);
          } else{
            dist = irotator->contact_distance(ior, isystem);
          }
          const std::vector<int>& idx = irotator->indices_[ior];
          contact.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], dist);
          //INFO(ior << " " << rot->contact_distance(ior, sys));
          if (!verbf.empty()) {
            #pragma omp critical
            {
              //const Position& pos = system->configuration().particle(1).site(0).position();
              pos.set_from_spherical({dist, irotator->stheta_[ior], irotator->sphi_[ior]});
              file << dist << "," << irotator->stheta_[ior] << "," << irotator->sphi_[ior]
                << "," << irotator->eulers_[ior].phi() << "," << irotator->eulers_[ior].theta()
                << "," << irotator->eulers_[ior].psi() << "," << pos.str() << std::endl;
            }
          }
        }
      }
    }
  }
  return contact;
}

bool BuildRecursiveTable::is_iso_(const bool is_2d, const int bsize) const {
  bool is_iso;
  if (is_2d) {
    is_iso = bsize == 1;
  } else {
    is_iso = bsize == 2;
  }
  DEBUG("is_iso: " << is_iso);
  return is_iso;
}

void BuildRecursiveTable::get_table_ptrs_(const bool is_cutoff, RecursiveTable * recur,
    RTable *& contact, RTable *& cutoff, RTable *& energy, RTable *& table) const {
  if (!recur) {
    INFO("detecting 1D iso-iso energy table");
    ASSERT(base_mrt_, "er");
    table = &base_mrt_->energy_[0][0];
  } else if (recur->energy_.size() > 0) {
    INFO("detecting 6D energy table");
    contact = &recur->contact_[0][0];
    cutoff = &recur->cutoff_[0][0];
    energy = &recur->energy_[0][0];
    table = energy;
  } else if (recur->cutoff_.size() > 0 && is_cutoff) {
    INFO("detecting 5D cutoff table");
    cutoff = &recur->cutoff_[0][0];
    table = cutoff;
  } else if (recur->contact_.size() > 0) {
    INFO("detecting 5D contact table");
    contact = &recur->contact_[0][0];
    table = contact;
  } else if (recur->energy3d_.size() > 0) {
    INFO("detecting 3D energy table");
    contact = &recur->contact2d_[0][0];
    cutoff = &recur->cutoff2d_[0][0];
    energy = &recur->energy3d_[0][0];
    table = energy;
  } else if (recur->cutoff2d_.size() > 0 && is_cutoff) {
    INFO("detecting 2D cutoff table");
    cutoff = &recur->cutoff2d_[0][0];
    table = cutoff;
  } else if (recur->contact2d_.size() > 0) {
    INFO("detecting 2D contact table");
    contact = &recur->contact2d_[0][0];
    table = contact;
  } else if (recur->energy2d_.size() > 0) {
    INFO("detecting 2D energy table");
    contact = &recur->contact1d_[0][0];
    cutoff = &recur->cutoff1d_[0][0];
    energy = &recur->energy2d_[0][0];
    table = energy;
  } else if (recur->cutoff1d_.size() > 0 && is_cutoff) {
    INFO("detecting 1D cutoff table");
    cutoff = &recur->cutoff1d_[0][0];
    table = cutoff;
  } else if (recur->contact1d_.size() > 0) {
    INFO("detecting 1D contact table");
    contact = &recur->contact1d_[0][0];
    table = contact;
  } else {
    FATAL("unrecognized data");
  }
  DEBUG("energy? " << energy);
}

void BuildRecursiveTable::base_iso_iso_(std::ostream& ostr, MonteCarlo * mc) {
  base_mrt_ = std::make_shared<ModelRecursiveTable>();
  DEBUG("build 1d table");
  double lower = NEAR_INFINITY, upper = cutoff_;
  { DEBUG("decide if upper is set by cutoff, or is to be determined.");
    bool find_upper = cutoff_ == -1;
    DEBUG("read training data line by line to find the global lower" <<
         "reading all lines required for every processor to use the same global bounds");
    std::ifstream mayer(mayer_training_file_);
    ASSERT(mayer.good(), "mayer_training_file: " << mayer_training_file_);
    bool lower_found = false;
    bool upper_found = false;
    std::string line;
    std::vector<std::string> sdat;
    while (std::getline(mayer, line)) {
      sdat = split(line, ',');
      const double dist = str_to_double(sdat[0]);
      const double en = str_to_double(sdat[1]);
      if (en < hard_limit_u_ && dist < lower) {
        lower = dist;
        lower_found = true;
      }
      if (find_upper && en != 0 && dist > upper) {
        upper = dist;
        upper_found = true;
      }
    }
    if (base_mrt_) {
      lower = std::sqrt(lower);
      if (find_upper) {
        upper = std::sqrt(upper);
      }
    }
    DEBUG("lower: " << lower << " upper: " << upper);
    ASSERT(lower_found, "lowest inner distance not found");
    ASSERT(!find_upper || upper_found, "upper distance not found");
  }

  RecursiveTable1D energy_table(argtype({{"num", str(num_z_)}}));
  build_table_({{lower, upper}}, &energy_table, mc);
  base_mrt_->lower_ = {{lower}};
  base_mrt_->upper_ = {{upper}};
  base_mrt_->energy_ = {{energy_table}};
  base_mrt_->serialize(ostr);
}

void BuildRecursiveTable::base(const bool is_energy, MonteCarlo * mc) {
  const Configuration& config = mc->system().configuration();
  std::stringstream ss;
  if (config.particle(0).num_sites() == 1) {
    base_iso_iso_(ss, mc);
  } else {
    const bool second_particle_iso = config.particle(1).num_sites() == 1;
    const bool is_2d  = static_cast<int>(config.dimension()) == 2;
    DEBUG("is_2d: " << is_2d);

    if (!is_energy) {
      base_rt_ = std::make_shared<RecursiveTable>();
    } else {
      if (!base_rt_) {
        ASSERT(!input_recursive_table_.empty(),
          "BuildRecursiveTable::input_recursive_table argument required");
        base_rt_ = std::make_shared<RecursiveTable>(RecursiveTable().from_file(input_recursive_table_));
      }
    }

    // prep systems for use with OMP
    int num_threads = get_num_threads_();
    INFO("num_threads " << num_threads);
    std::vector<System> systems(num_threads - 1);
    for (int thread = 0; thread < num_threads - 1; ++thread) {
      systems[thread] = deep_copy(mc->system());
    }

    std::string verbf = verbose_name_(0);
    if (is_energy) {
      verbf += "en";
    }
    const std::vector<std::vector<double> >& abnd = Rotator().gen_global_bounds(config);
    if (second_particle_iso && is_2d) {
      if (!is_energy) {
        const auto contact = build_1dcontact_(abnd, verbf, mc->get_system(), &systems);
        base_rt_->contact1d_ = {{contact}};
        if (!contact_only_) {
          const auto cutoff = build_1dcontact_(abnd, verbf+"cut", mc->get_system(), &systems, true);
          base_rt_->cutoff1d_ = {{cutoff}};
        }
      } else {
        ASSERT(!contact_only_, "err");
        const auto& contact = base_rt_->contact1d_[0][0];
        const auto& cutoff  = base_rt_->cutoff1d_[0][0];
        const auto energy = build_2denergy_(abnd, {0., 1.}, verbf, contact, cutoff, mc->get_system(), &systems);
        base_rt_->energy2d_ = {{energy}};
      }
    } else if (( second_particle_iso && !is_2d) ||
               (!second_particle_iso &&  is_2d)) {
      if (!is_energy) {
        const auto contact = build_2dcontact_(abnd, verbf, mc->get_system(), &systems);
        base_rt_->contact2d_ = {{contact}};
        if (!contact_only_) {
          const auto cutoff = build_2dcontact_(abnd, verbf+"cut", mc->get_system(), &systems, true);
          base_rt_->cutoff2d_ = {{cutoff}};
        }
      } else {
        ASSERT(!contact_only_, "err");
        ASSERT(base_rt_, "error");
        ASSERT(base_rt_->contact2d_.size() > 0, "error");
        const auto& contact = base_rt_->contact2d_[0][0];
        ASSERT(base_rt_->cutoff2d_.size() > 0, "error");
        const auto& cutoff  = base_rt_->cutoff2d_[0][0];
        const auto energy = build_3denergy_(abnd, {0., 1.}, verbf, contact, cutoff, mc->get_system(), &systems);
        base_rt_->energy3d_ = {{energy}};
      }
    } else {
      if (!is_energy) {
        const auto contact = build_contact_(abnd, verbf, mc->get_system(), &systems);
        base_rt_->contact_ = {{contact}};
        if (!contact_only_) {
          INFO("building cutoff");
          const auto cutoff = build_contact_(abnd, verbf+"cut", mc->get_system(), &systems, true);
          base_rt_->cutoff_ = {{cutoff}};
        }
      } else {
        ASSERT(!contact_only_, "err");
        ASSERT(base_rt_->contact_.size() > 0, "error");
        const auto& contact = base_rt_->contact_[0][0];
        ASSERT(base_rt_->cutoff_.size() > 0, "error");
        const auto& cutoff  = base_rt_->cutoff_[0][0];
        const auto energy = build_energy_(abnd, {0., 1.}, verbf, contact, cutoff, mc->get_system(), &systems);
        base_rt_->energy_ = {{energy}};
      }
    }
    base_rt_->serialize(ss);
  }

  // write serialization to file
  std::string outfilename = output_file_;
  if (stage_ == "all") {
    outfilename += "_base";
  }
  ASSERT(!outfilename.empty(), "Error");
  std::ofstream outfile(outfilename);
  ASSERT(outfile.good(), "Error");
  outfile << ss.str();
}

void BuildRecursiveTable::train(const bool is_cutoff, const bool is_energy, MonteCarlo * mc) {
  INFO("Begin training");
  ASSERT(!is_cutoff || !contact_only_, "err");
  if (!base_table_file_.empty()) {
    base_rt_ = std::make_shared<RecursiveTable>(RecursiveTable().from_file(base_table_file_));
  }
  RTable * table, * contact = NULL, * cutoff = NULL, * energy = NULL;
  get_table_ptrs_(is_cutoff, base_rt_.get(), contact, cutoff, energy, table);
  INFO("table dim " << table->dimension());
  int read_en_index = table->dimension() + 1;
  if (is_energy) {
    --read_en_index;
  }
  INFO("read_en_index " << read_en_index);
  //ASSERT(energy == NULL, "if energy, needs to train nested for contact and cutoff separately first, then use those nested to train energy");

  // initialize extra verbose files
  const Configuration& config = mc->configuration();
  const std::vector<std::vector<double> > abnd = Rotator().gen_global_bounds(config);
  Position pos, pos2;
  bool is_2d = false, is_iso = false;
  std::ofstream file, file2;
  std::string verbf = verbose_name_(1);
  if (is_energy) {
    verbf += "en";
  } else if (is_cutoff) {
    verbf += "cut";
  }
  if (base_mrt_) {
    if (!verbf.empty()) {
      file.open(verbf);
    }
  } else {
    is_2d = config.dimension() == 2;
    INFO("is_2d? " << is_2d);
    is_iso = config.particle(1).num_sites() == 1;
    INFO("is_iso? " << is_iso);
    if (extra_verbose_) {
      file.open(verbf+"an");
      file2.open(verbf+"vis");
      if (is_2d) {
        if (is_iso) {
          file << "x,y" << std::endl;
        } else {
          file << "x,y,phi" << std::endl;
        }
      } else {
        file << "x,y,z" << std::endl;
      }
    }
  }

  // read mayer training file line by line until eof
  // parallelize by chunks of lines and not every other line,
  // because sampling results in correlated data, which could lead to
  // the same max_bins in every processor.
  std::ifstream mayer(mayer_training_file_);
  ASSERT(mayer.good(), "mayer_training_file: " << mayer_training_file_);
  std::string line;
  std::vector<std::string> sdat;
  std::vector<double> scaled_coords(table->dimension());
  std::vector<double> scaled_orient(table->dimension()-1);
  rbin max_bins;
  max_bins.first.resize(table->dimension());
  rbins_ = std::make_shared<std::vector<rbin> >();
  const double numl = num_lines(mayer_training_file_);
  int iline = 0;
  bool reading = true;
  while (std::getline(mayer, line) && reading) {
    int current_proc = int(num_processors_*iline/numl);
    ASSERT(current_proc >= 0 && current_proc < num_processors_,
      "current_proc: " << current_proc << "num_proc:" << num_processors_);
    if (current_proc > processor_) {
      reading = false;
    } else if (current_proc == processor_) {
      sdat = split(line, ',');
      //INFO("sdat " << feasst_str(sdat));
      // WARN assumes data is r, {orientations}, energy
      ASSERT(((contact || cutoff) && static_cast<int>(sdat.size()) == table->dimension() + 2) ||
             (energy && static_cast<int>(sdat.size()) == table->dimension() + 1) ||
             (table->dimension() == 1 && static_cast<int>(sdat.size() == 2)),
             "data not recognized: " << line);
      double dist = str_to_double(sdat[0]);
      DEBUG("dist " << dist);

      // obtain en and scaled coordinates
      double en, z=-1;
      if (base_mrt_) {
        ASSERT(!energy, "Er");
        ASSERT(!contact, "Er");
        ASSERT(!cutoff, "Er");
        en = str_to_double(sdat[1]);
        dist = std::sqrt(dist);
        DEBUG("updated dist " << dist);
        const double rh = base_mrt_->lower_[0][0];
        const double rc = base_mrt_->upper_[0][0];
        z = (dist - rh)/(rc - rh);
        if (z >= 0. && z <= 1.) {
          scaled_coords[0] = z;
        } else {
          ++iline;
          continue;
        }
      } else {
        en = str_to_double(sdat[read_en_index]);
        if (energy) {
          for (int ist = 0; ist < table->dimension() - 1; ++ist) {
            scaled_orient[ist] = str_to_double(sdat[ist+1]);
          }
          const double rh = contact->linear_interpolation(scaled_orient);
          const double rc = cutoff->linear_interpolation(scaled_orient);
          const double z = (dist - rh)/(rc - rh);
          if (z >= 0. && z <= 1.) {
            for (int ist = 0; ist < table->dimension() - 1; ++ist) {
              scaled_coords[ist] = scaled_orient[ist];
            }
            scaled_coords[table->dimension() - 1] = z;
          } else {
            ++iline;
            continue;
          }
        } else {
          for (int ist = 0; ist < table->dimension(); ++ist) {
            scaled_coords[ist] = str_to_double(sdat[ist+1]);
          }
        }
      }
      DEBUG("en " << en);
      DEBUG("scaled_coords: " << feasst_str(scaled_coords));
      const double interp = table->linear_interpolation(scaled_coords);
      DEBUG("interp: " << interp);

      // obtain criteria
      double criteria = 0;
      bool too_big = false, too_small = false;
      if (energy) {
        DEBUG("obtaining energy criteria");
        criteria = std::abs(std::exp(-beta_*en) - std::exp(-beta_*interp));
      } else if (is_cutoff) {
        DEBUG("obtaining cutoff criteria");
        too_small = en != 0. && dist > interp;
        too_big   = en == 0. && dist < interp;
        if (too_small || too_big) {
          criteria = std::abs(dist - interp);
        }
      } else if (contact) {
        DEBUG("obtaining contact criteria");
        too_small = en > hard_limit_u_ && dist > interp;
        too_big   = en < hard_limit_u_ && dist < interp;
        if (too_small || too_big) {
          criteria = std::abs(dist - interp);
        }
      } else {
        DEBUG("obtaining iso-iso 1d criteria");
        ASSERT(base_mrt_, "Assumed ModelRecursiveTable");
        criteria = dist*dist*std::abs(std::exp(-beta_*en) - std::exp(-beta_*interp));
      }
      DEBUG("criteria " << criteria << " min_crit " << min_criteria_);

      // print verbose file
      if (file.is_open()) {//!verbf.empty() && extra_verbose_) {
        if (base_mrt_) {
          file << feasst_str(std::vector<double>({dist, z, en, interp, en-interp, criteria})) << std::endl;
//          file << dist << "," << z << "," << en << "," << interp << "," << en - interp << "," << criteria << std::endl;
        } else {
          if (is_2d) {
            const double dtht = abnd[0][1]-abnd[0][0];
            if (is_iso) {
              pos.set_from_spherical({interp, dtht*scaled_coords[0]});
              file << pos.str() << std::endl;
              pos2.set_from_spherical({dist, dtht*scaled_coords[0]});
              file2 << pos2.str() << "," << en << "," << too_small << "," << too_big << std::endl;
            } else {
              pos.set_from_spherical({interp, dtht*scaled_coords[0]});
              file << pos.str() << "," << (abnd[1][1] - abnd[1][0])*(scaled_coords[1] - 0.5) << std::endl;
            }
          } else {
            pos.set_from_spherical({interp, (abnd[0][1]-abnd[0][0])*scaled_coords[0],
                                            (abnd[1][1]-abnd[1][0])*scaled_coords[1]});
            file << "0 " << pos.coord(0) << " " << pos.coord(1) << " " << pos.coord(2) << std::endl;
          }
        }
      }

      if ( (!is_energy && criteria > min_criteria_) ||
           ( is_energy && criteria > min_criteria_energy_) ) {
        DEBUG("criteria " << criteria << " > min: " << min_criteria_);
        for (int rdim = 0; rdim < static_cast<int>(max_bins.first.size()); ++rdim) {
          max_bins.first[rdim] = table->value_to_lowest_bin(rdim, scaled_coords[rdim]);
        }
        max_bins.second = criteria;
        DEBUG("max_bins " << feasst_str(max_bins));
        int index;
        if (find_in_list(max_bins.first, *rbins_, &index)) {
        //if (find_in_list_of_pair(max_bins, *rbins_, &index)) {
          DEBUG("Found, updating max crit");
          if (criteria > (*rbins_)[index].second) {
            (*rbins_)[index].second = criteria;
          }
        } else {
          DEBUG("Not found. Adding to rbins");
          rbins_->push_back(max_bins);
        }
      }
    }
    ++iline;
  }
  INFO("num nested " << rbins_->size());
  INFO("nested percentage: " << static_cast<double>(rbins_->size())/table->num());

  std::string outfilename = output_file_;
  if (stage_ == "all") {
    outfilename += "_trained";
  }
  if (is_energy) {
    outfilename += "_en";
  } else if (is_cutoff) {
    outfilename += "_cut";
  }
  sort_and_print(rbins_.get(), outfilename);
}

void sort_and_print(std::vector<rbin> * rbins, const std::string& filename) {
  DEBUG("sort by criteria");
  std::sort(rbins->begin(), rbins->end(), [](const rbin& a, const rbin& b) {
    return a.second > b.second;
  });
  std::ofstream outfile(filename);
  ASSERT(outfile.good(), "cannot initialize: " << filename);
  outfile << "#[bins,]:criteria" << std::endl;
  for (const rbin& rb : *rbins) {
    outfile << feasst_str(rb.first) << ":" << MAX_PRECISION << rb.second << std::endl;
  }
}

std::vector<rbin> read_rbins(const std::string& filename, const int processor, const int num_processors) {
  DEBUG("get rbins from trained file");
  std::vector<rbin> rbins;
  const int numl = num_lines(filename);
  ASSERT(numl >= 0, "numl: " << numl);
  std::ifstream trained(filename);
  ASSERT(trained.good(), "Error reading trained_file: " << filename);
  rbin bin;
  rbins.resize(numl);
  std::string line;
  std::vector<std::string> sdat1, sdat2;
  int iline = 0;
  int inum = 0;
  while (std::getline(trained, line)) {
    if (line.front() != '#' && !line.empty()) {
      sdat1 = split(line, ':');
      ASSERT(static_cast<int>(sdat1.size()) == 2, sdat1.size());
      sdat2 = split(sdat1[0], ',');
      bin.first.resize(static_cast<int>(sdat2.size()));
      for (int ist = 0; ist < static_cast<int>(sdat2.size()); ++ist) {
        bin.first[ist] = str_to_int(sdat2[ist]);
      }
      bin.second = str_to_double(sdat1[1]);
      ASSERT(inum < static_cast<int>(rbins.size()),
        "inum: " << inum << " >= rbins.size:" << rbins.size());
      rbins[inum] = bin;
      ++inum;
    }
    ++iline;
    ASSERT(iline <= numl, "error");
  }
  rbins.resize(inum);
  return rbins;
}

void BuildRecursiveTable::build(const bool is_cutoff, const bool is_energy, MonteCarlo * mc) {
  INFO("building is_cutoff? " << is_cutoff << " is_energy? " << is_energy);
  ASSERT(!is_cutoff || !contact_only_, "err");
  if (!rbins_) {
    rbins_ = std::make_shared<std::vector<rbin> >();
    ASSERT(!trained_file_.empty(), "BuildRecursiveTable::trained_file is required.");
    *rbins_ = read_rbins(trained_file_, processor_, num_processors_);
  }
  if (rbins_->size() == 0) {
    INFO("Nothing found without criteria, so no recursions.");
  } else {

    // global rotation bounds
    const Configuration& config = mc->system().configuration();
    std::vector<std::vector<double> > abnd;
    if (base_mrt_) {
      abnd = {{base_mrt_->lower_[0][0], base_mrt_->upper_[0][0]}};
    } else {
      abnd = Rotator().gen_global_bounds(config);
    }
    DEBUG("abnd: " << feasst_str(abnd));
    std::vector<double> zbnds(2);

    // tables
    RecursiveTable * rt = NULL;
    if (!base_mrt_ && !base_rt_) {
      base_rt_ = std::make_shared<RecursiveTable>(RecursiveTable().from_file(base_table_file_));
      rt = base_rt_.get();
    } else if (base_rt_) {
      rt = base_rt_.get();
    }

    RTable * contact, * table;
    RTable * cutoff = NULL;
    RTable * energy = NULL;
    get_table_ptrs_(is_cutoff, rt, contact, cutoff, energy, table);

    // nested table
    RecursiveTable1D nested1d;
    RecursiveTable2D nested2d;
    RecursiveTable3D nested3d;
    RecursiveTable5D nested5d;
    RecursiveTable6D nested6d;
    RTable * nested = NULL;
    DEBUG("table->dimension() " << table->dimension());
    if (base_mrt_) {
      nested1d = RecursiveTable1D(argtype({{"num", str(num_z_)}}));
    }
    if (table->dimension() == 1) {
      nested = &nested1d;
    } else if (table->dimension() == 2) {
      nested = &nested2d;
    } else if (table->dimension() == 3) {
      nested = &nested3d;
    } else if (table->dimension() == 5) {
      nested = &nested5d;
    } else if (table->dimension() == 6) {
      nested = &nested6d;
    } else {
      FATAL("table");
    }

    DEBUG("rbins " << feasst_str(*rbins_));
    std::vector<std::vector<double> > bounds;
    if (is_energy) {
      resize((*rbins_)[0].first.size() - 1, 2, &bounds);
    } else {
      resize((*rbins_)[0].first.size(), 2, &bounds);
    }

    // copy systems for use with parallelization
    std::vector<System> systems;
    { int num_threads = get_num_threads_();
      DEBUG("num_threads " << num_threads);
      systems.resize(num_threads - 1);
      for (int thread = 0; thread < num_threads - 1; ++thread) {
        systems[thread] = deep_copy(mc->system());
      }
    }

    std::string verbf;
    for (int ibin = 0; ibin < static_cast<int>(rbins_->size()); ++ibin) {
      DEBUG("ibin " << ibin);
      //for (const std::vector<int>& bins : rbins) {
      const std::vector<int>& bins = (*rbins_)[ibin].first;
      DEBUG("bins: " << feasst_str(bins));
      verbf = verbose_name_(1 + ibin*num_processors_ + processor_);
      if (is_energy) {
        verbf += "en";
      } else if (is_cutoff) {
        verbf += "cut";
      }
      int num_rdim = static_cast<int>(bins.size());
      if (is_energy) {
        --num_rdim;
        zbnds[0] = table->bin_to_value(num_rdim, bins[num_rdim]);
        zbnds[1] = table->bin_to_value(num_rdim, bins[num_rdim] + 1);
        DEBUG("zbnds: " << feasst_str(zbnds));
      }
      for (int rdim = 0; rdim < num_rdim; ++rdim) {
        const std::vector<double>& bn = abnd[rdim];
        const double lr = bn[0] + (bn[1] - bn[0])*table->bin_to_value(rdim, bins[rdim]);
        const double ur = bn[0] + (bn[1] - bn[0])*table->bin_to_value(rdim, bins[rdim] + 1);
        bounds[rdim][0] = lr;
        bounds[rdim][1] = ur;
      //DEBUG("rdim " << rdim << " bn " << feasst_str(bn) << " table " << table->bin_to_value(rdim, bins[rdim]) << " table+1 " << table->bin_to_value(rdim, bins[rdim]+1));
      }
      DEBUG("bounds: " << feasst_str(bounds));
      if (base_mrt_) {
        build_table_(bounds, &nested1d, mc);
      } else if (table->dimension() == 1) {
        nested1d = build_1dcontact_(bounds, verbf, mc->get_system(), &systems, is_cutoff);
      } else if (table->dimension() == 2) {
        if (is_energy) {
          ASSERT(!is_cutoff, "err");
          const RecursiveTable1D& cntct = base_rt_->contact1d_[0][0];
          const RecursiveTable1D& cutff = base_rt_->cutoff1d_[0][0];
          nested2d = build_2denergy_(bounds, zbnds, verbf, cntct, cutff, mc->get_system(), &systems);
        } else {
          nested2d = build_2dcontact_(bounds, verbf, mc->get_system(), &systems, is_cutoff);
        }
      } else if (table->dimension() == 3) {
        ASSERT(is_energy && !is_cutoff, "err");
        const RecursiveTable2D& cntct = base_rt_->contact2d_[0][0];
        const RecursiveTable2D& cutff = base_rt_->cutoff2d_[0][0];
        nested3d = build_3denergy_(bounds, zbnds, verbf, cntct, cutff, mc->get_system(), &systems);
      } else if (table->dimension() == 5) {
        ASSERT(!is_energy, "err");
        nested5d = build_contact_(bounds, verbf, mc->get_system(), &systems, is_cutoff);
      } else if (table->dimension() == 6) {
        ASSERT(is_energy && !is_cutoff, "err");
        const RecursiveTable5D& cntct = base_rt_->contact_[0][0];
        const RecursiveTable5D& cutff = base_rt_->cutoff_[0][0];
        nested6d = build_energy_(bounds, zbnds, verbf, cntct, cutff, mc->get_system(), &systems);
      }
      DEBUG("nested size " << nested->dimension() << " num " << nested->num());
      ASSERT(nested->num() > 1, "empty nested");
      ASSERT(static_cast<int>(bins.size()) > 0, "bin size: " << bins.size());
      ASSERT(static_cast<int>(bins.size()) == table->dimension(), "bin size: " << bins.size());
      table->insert(bins, *nested);
    }
  }
  std::stringstream ss;
  if (base_rt_) {
    base_rt_->serialize(ss);
  } else if (base_mrt_) {
    base_mrt_->serialize(ss);
  }
  ASSERT(!output_file_.empty(), "Error");
  std::ofstream outfile(output_file_);
  ASSERT(outfile.good(), "Error");
  outfile << ss.str();
}

void BuildRecursiveTable::run(MonteCarlo * mc) {
  INFO("stage " << stage_);
  if (stage_ == "base") {
    base(energy_only_, mc);
  } else if (stage_ == "train") {
    if (energy_only_) {
      train(false, true, mc);
    } else {
      train(false, false, mc);
      if (!contact_only_) {
        train(true, false, mc);
      }
    }
  } else if (stage_ == "build_contact") {
      build(false, false, mc);
  } else if (stage_ == "build_cutoff") {
      build(true, false, mc);
  } else if (stage_ == "build_energy") {
      build(false, true, mc);
//    if (energy_only_) {
//      build(false, true, mc);
//    } else {
//      build(false, false, mc);
//      if (!contact_only_) {
//        build(true, false, mc);
//      }
//    }
  } else if (stage_ == "all") {
    base(false, mc);
    train(false, false, mc);
    build(false, false, mc);
    if (!contact_only_) {

      // build cutoff
      train(true, false, mc);
      build(true, false, mc);

      // build energy
      base(true, mc);
      train(false, true, mc);
      build(false, true, mc);
    }
  } else {
    FATAL("unrecognized stage: " << stage_);
  }
}

}  // namespace feasst
