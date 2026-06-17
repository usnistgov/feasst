#include <fstream>
#include "threads/include/thread_omp.h"
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
#include "aniso/include/model_recursive_table.h"
#include "aniso/include/recursive_table.h"
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
  feasst_deserialize(&extra_verbose_, istr);
}

void BuildRecursiveTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
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
  feasst_serialize(extra_verbose_, ostr);
}

void BuildRecursiveTable::build_table_(const std::vector<std::vector<double> >& bounds, const std::vector<std::vector<double> >& data, RecursiveTable1D * tab, MonteCarlo * mc) const {
  DEBUG("building table. inner:" << bounds[0][0] << " outer:" << bounds[0][1]);
  DEBUG("data size:" << data.size());
  DEBUG("data size:" << data[0].size());
  if (data[0].size() == 2) {
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
      const double z = static_cast<double>(i)/(tab->num() - 1);
      const double r = lower + z*(upper - lower);
      DEBUG("r:" << r);
      coords[1][0] = r;
      sys->get_configuration()->update_positions(coords);
      const double en = sys->energy();
      tab->set_data(i, en);
    }
    sys->get_configuration()->update_positions(old_coords);
    sys->energy();
  } else {
    FATAL("unrecognized data size:" << data[0].size());
  }
  DEBUG("table built");
}

int BuildRecursiveTable::analyze_table_(const double lower, const double upper, const std::vector<std::vector<double> >& data, const Table1D& table, const double beta, const std::string& filename, double * max_criteria) {
  // compare the table with the data
  // find the global max of |e^-bU - e^-bU_table| to insert a recurive table at that point
  double max_z = -1;
  *max_criteria = -1;
  std::ofstream file(filename);
  const double zfac = 1./(upper - lower);
  for (int idat = 0; idat < static_cast<int>(data.size()); ++idat) {
    const std::vector<double>& dat = data[idat];
    const double r = std::sqrt(dat[0]);
    if (r >= lower && r <= upper) {
      const double en = dat[1];
      double z = (r - lower)*zfac;
      DEBUG("z " << MAX_PRECISION << z);
      if (z <= 0.) {
        ASSERT(z + NEAR_ZERO >= 0, "z:" << MAX_PRECISION << z);
        z = NEAR_ZERO;
      }
      if (z >= 1.) {
        ASSERT(z - NEAR_ZERO <= 1, "z:" << MAX_PRECISION << z);
        z = 1. - NEAR_ZERO;
      }
      DEBUG("z " << MAX_PRECISION << z);
      const double entab = table.linear_interpolation(z);
      const double crit = r*r*std::abs(std::exp(-beta*en) - std::exp(-beta*entab));
      if (!filename.empty()) {
        file << r << " " << z << " " << en << " " << entab << " " << en - entab << " " << crit << std::endl;
      }
      if (crit > *max_criteria) {
        max_z = z;
        *max_criteria = crit;
      }
      //INFO("r2: " << r2 << " z: " << z << " en: " << en << " tab: " << entab << " diff: " << en - entab);
    }
  }
  ASSERT(*max_criteria > 0, "max_criteria: " << *max_criteria);
  INFO("max_z:" << max_z << " max_criteria:" << *max_criteria);
  const int max_bin = table.value_lowest_bin(max_z);
  INFO("max_bin:" << max_bin);
  return max_bin;
}


std::string BuildRecursiveTable::verbose_name_(const int iteration) const {
  if (verbose_file_.empty()) {
    return std::string("");
  }
  return verbose_file_ + "_" + str(iteration);
}

void read_data_(const std::string& filename, std::vector<std::vector<double> > * data) {
  DEBUG("obtain data from file");
  data->clear();
  ASSERT(!filename.empty(), "requires mayer_training_file");
  std::ifstream file(filename);
  ASSERT(file.good(), "Cannot find file: " << filename);
  std::string line;
  while (!file.eof()) {
    std::getline(file, line);
    if (!line.empty()) {
      std::vector<std::string> dat = split(line, ',');
      std::vector<double> datd;
      for (const std::string& st : dat) {
        const double val = str_to_double(st);
        datd.push_back(val);
      }
      if (datd.size() > 0) {
        data->push_back(datd);
      }
    }
  }
  DEBUG("sorting");
  std::sort(data->begin(), data->end());
  DEBUG("sorted");
  DEBUG("remove duplicants / redundant");
  auto it = std::unique(data->begin(), data->end());
  data->erase(it, data->end());

//  // print data to screen
//  DEBUG(data.size());
//  for (const auto& d : data) {
//    DEBUG(feasst_str(d));
//  }
  file.close();
}

double global_lower_(const double hard_limit_u, std::vector<std::vector<double> > * data) {
  DEBUG("find lowest distance with energy less than hard_limit. Remove data below");
  int found = -1;
  double lower = -1;
  for (int index = 0; index < static_cast<int>(data->size()); ++index) {
    const std::vector<double>& en = (*data)[index];
    if (en[1] < hard_limit_u) {
      lower = std::sqrt(en[0]);
      found = index;
      DEBUG("found inner distance of " << en[0]);
      break;
    }
  }
  ASSERT(found != -1, "lowest inner distance not found");
  data->erase(data->begin(), data->begin() + found);
  return lower;
}

double global_upper_(const double cutoff, std::vector<std::vector<double> > * data) {
  double upper = cutoff;
  if (upper == -1) {
    DEBUG("find upper distance with energy equal to zero.");
    int found = -1;
    for (int index = static_cast<int>(data->size()) - 1; index >= 0; --index) {
      const std::vector<double>& en = (*data)[index];
      if (en[1] != 0) {
        found = index;
        DEBUG("found outer distance of " << en[0]);
        break;
      }
    }
    ASSERT(found != -1, "largest outer distance not found");
    DEBUG("found : " << found << " max " << data->back()[0]);
    DEBUG("found+1 : " << (*data)[found][0] << " " << (*data)[found+1][0]);
    if (found < static_cast<int>(data->size()) - 1) {
      upper = std::sqrt((*data)[found + 1][0]);
      data->erase(data->begin() + found + 1, data->end());
    } else {
      upper = std::sqrt(data->back()[0]);
    }
    DEBUG("upper: " << MAX_PRECISION << upper);
  } else {
    DEBUG("remove data above cutoff");
    int found = -1;
    for (int index = static_cast<int>(data->size()) - 1; index >= 0; --index) {
      const std::vector<double>& en = (*data)[index];
      if (en[0]*en[0] < cutoff) {
        found = index;
        break;
      }
    }
    if (found < static_cast<int>(data->size()) - 1) {
      data->erase(data->begin() + found + 1, data->end());
    }
  }
  return upper;
}

RecursiveTable1D BuildRecursiveTable::build_1dcontact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system, const bool cutoff) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
  RecursiveTable1D contact(argtype({{"num", str(rotator.sizes_[0])}}));
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf);
    file << "r,theta,x,y" << std::endl;
  }
  Position pos(argtype({{"dimension", "2"}}));
  INFO("num " << rotator.num_orientations());
  for (int ior = 0; ior < rotator.num_orientations(); ++ior) {
    double dist;
    if (cutoff) {
      dist = rotator.cutoff_distance(ior, system);
    } else{
      dist = rotator.contact_distance(ior, system);
    }
    const std::vector<int>& idx = rotator.indices_[ior];
    contact.set_data(idx[0], dist);
    if (!verbf.empty()) {
      pos.set_from_spherical({dist, rotator.stheta_[ior]});
      file << dist << "," << rotator.stheta_[ior] << "," << pos.str()
           << std::endl;
    }
  }
  return contact;
}

RecursiveTable2D BuildRecursiveTable::build_2dcontact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system, const bool cutoff) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
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
  Position pos(argtype({{"dimension", str(dimen)}}));
  INFO("num " << rotator.num_orientations());
  for (int ior = 0; ior < rotator.num_orientations(); ++ior) {
    double dist;
    if (cutoff) {
      dist = rotator.cutoff_distance(ior, system);
    } else{
      dist = rotator.contact_distance(ior, system);
    }
    DEBUG("dist " << dist);
    const std::vector<int>& idx = rotator.indices_[ior];
    DEBUG("idx " << feasst_str(idx));
//    // HWH DEBUG TESTING FOLLOWS, CAN DELETE
//    if (std::abs(rotator.eulers_[ior].phi() - PI) < 1e-6) {
//      INFO("theta " << rotator.stheta_[ior] << " phi " << rotator.eulers_[ior].phi() << " idx " << idx[1] << " dist " << dist);
//    }
    contact.set_data(idx[0], idx[1], dist);
    if (!verbf.empty()) {
      if (dimen == 3) {
        pos.set_from_spherical({dist, rotator.stheta_[ior], rotator.sphi_[ior]});
        file << dist << "," << rotator.stheta_[ior] << "," << rotator.sphi_[ior]
          << "," << pos.str() << std::endl;
      } else {
        pos.set_from_spherical({dist, rotator.stheta_[ior]});
        file << dist << "," << rotator.stheta_[ior] << "," << rotator.eulers_[ior].phi()
          << "," << pos.str() << std::endl;
      }
    }
  }
  return contact;
}

// HWH copied from below
RecursiveTable6D BuildRecursiveTable::build_energy_(const std::vector<std::vector<double> >& bounds, const std::vector<double>& zbnds, const std::string& verbf, const RecursiveTable5D& contact, const RecursiveTable5D& cutoff, System * system) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
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
  Position pos(argtype({{"dimension", str(dimen)}}));
  INFO("num " << rotator.num_orientations()*num_z_);
  double s1, s2, s3, s4, s5;
  const double dzpn = (zbnds[1] - zbnds[0])/static_cast<double>(num_z_ - 1);
  bool full_range = false;
  int subsize = num_z_;
  if (std::abs(zbnds[1] - zbnds[0] - 1) < 1e-8) {
    full_range = true;
    --subsize;
  }
  INFO("full range " << full_range << " subsize " << subsize << " dzpn " << dzpn);
  for (int ior = 0; ior < rotator.num_orientations(); ++ior) {
    const std::vector<int>& idx = rotator.indices_[ior];
    const Euler& eul = rotator.eulers_[ior];
    scaled_relative_orientation(rotator.stheta_[ior], rotator.sphi_[ior], eul.phi(), eul.theta(), eul.psi(), &s1, &s2, &s3, &s4, &s5);
    //INFO("s1 " << s1 << " s2 " << s2);
    const double rh = contact.linear_interpolation(s1, s2, s3, s4, s5);
    const double rc = cutoff.linear_interpolation(s1, s2, s3, s4, s5);
    //INFO("rh " << rh);
    //INFO("rc " << rc);
    for (int idis = 0; idis < subsize; ++idis) {
      const double z = zbnds[0] + dzpn*static_cast<double>(idis);
      const double dist = rh + z*(rc - rh);
      double ener = rotator.energy(ior, dist, system);
      if (ener > NEAR_INFINITY_FLOAT) {
        ener = NEAR_INFINITY_FLOAT;
      }
      en.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], idis, ener);
      if (!verbf.empty()) {
        pos.set_from_spherical({dist, rotator.stheta_[ior], rotator.sphi_[ior]});
        file << dist << "," << rotator.stheta_[ior] << "," << rotator.sphi_[ior]
          << "," << eul.phi() << "," << eul.theta() << "," << eul.psi() << ","
          << pos.str() << "," << ener << std::endl;
      }
    }
    if (full_range) {
      // Set energy to zero at cutoff (z=1)
      en.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], num_z_ - 1, 0.);
    }
  }
  return en;
}

RecursiveTable3D BuildRecursiveTable::build_3denergy_(const std::vector<std::vector<double> >& bounds, const std::vector<double>& zbnds, const std::string& verbf, const RecursiveTable2D& contact, const RecursiveTable2D& cutoff, System * system) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
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
  Position pos(argtype({{"dimension", str(dimen)}}));
  INFO("num " << rotator.num_orientations()*num_z_);
  double s1, s2;
  const double dzpn = (zbnds[1] - zbnds[0])/static_cast<double>(num_z_ - 1);
  bool full_range = false;
  int subsize = num_z_;
  if (std::abs(zbnds[1] - zbnds[0] - 1) < 1e-8) {
    full_range = true;
    --subsize;
  }
  INFO("full range " << full_range << " subsize " << subsize << " dzpn " << dzpn);
  for (int ior = 0; ior < rotator.num_orientations(); ++ior) {
    const std::vector<int>& idx = rotator.indices_[ior];
    if (dimen == 3) {
      scaled_relative_orientation(rotator.stheta_[ior], rotator.sphi_[ior], dimen, &s1, &s2);
    } else if (dimen == 2) {
      scaled_relative_orientation(rotator.stheta_[ior], rotator.eulers_[ior].phi(), dimen, &s1, &s2);
    } else {
      FATAL("unrecognized dimen");
    }
    //INFO("s1 " << s1 << " s2 " << s2);
    // HWH I'm not sure why, but without this PBC the high end of the boundary doesn't use the recursive tabled rc and rh to determine the energy distances
//    double rh, rc;
//    if (std::abs(s2 - 1.) < 1e-6) {
//      rh = contact.linear_interpolation(s1, 0.);
//      rc = cutoff.linear_interpolation(s1, 0.);
//    } else {
//      rh = contact.linear_interpolation(s1, s2);
//      rc = cutoff.linear_interpolation(s1, s2);
//    }
    const double rh = contact.linear_interpolation(s1, s2);
    const double rc = cutoff.linear_interpolation(s1, s2);
    //INFO("rh " << rh);
    //INFO("rc " << rc);
    for (int idis = 0; idis < subsize; ++idis) {
      const double z = zbnds[0] + dzpn*static_cast<double>(idis);
      const double dist = rh + z*(rc - rh);
      double ener = rotator.energy(ior, dist, system);
      if (ener > NEAR_INFINITY_FLOAT) {
        ener = NEAR_INFINITY_FLOAT;
      }
//      // HWH DEBUG TESTING FOLLOWS, CAN DELETE
//      if (std::abs(rotator.eulers_[ior].phi() - PI) < 1e-6) {
//      //if (std::abs(rotator.eulers_[ior].phi() - PI) < 1e-6 ||
//      //    std::abs(rotator.eulers_[ior].phi() + PI) < 1e-6) {
//        INFO("theta " << rotator.stheta_[ior] << " s1 " << s1 << " phi " << rotator.eulers_[ior].phi() << " s2 " << s2 << " idx " << idx[1] << " rh " << rh << " rhflip " << contact.linear_interpolation(s1, 1.-s2) << " rc " << rc << " rcflip " << cutoff.linear_interpolation(s1, 1.-s2));
//      }
      en.set_data(idx[0], idx[1], idis, ener);
      if (!verbf.empty()) {
        file << dist << "," << rotator.stheta_[ior];
        if (dimen == 3) {
          pos.set_from_spherical({dist, rotator.stheta_[ior], rotator.sphi_[ior]});
          file << "," << rotator.sphi_[ior];
        } else if (dimen == 2) {
          pos.set_from_spherical({dist, rotator.stheta_[ior]});
          file << "," << rotator.eulers_[ior].phi();
        } else {
          FATAL("unrecognized dimen");
        }
        file << "," << pos.str() << "," << ener << std::endl;
      }
    }
    if (full_range) {
      // Set energy to zero at cutoff (z=1)
      en.set_data(idx[0], idx[1], num_z_ - 1, 0.);
    }
  }
  return en;
}

// HWH copied from above
RecursiveTable2D BuildRecursiveTable::build_2denergy_(const std::vector<std::vector<double> >& bounds, const std::vector<double>& zbnds, const std::string& verbf, const RecursiveTable1D& contact, const RecursiveTable1D& cutoff, System * system) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
  RecursiveTable2D en(argtype({{"num0", str(rotator.sizes_[0])},
                               {"num1", str(num_z_)}}));
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf);
    file << "r,theta,x,y,en" << std::endl;
  }
  Position pos(argtype({{"dimension", "2"}}));
  INFO("num " << rotator.num_orientations()*num_z_);
  double s1;
  const double dzpn = (zbnds[1] - zbnds[0])/static_cast<double>(num_z_ - 1);
  bool full_range = false;
  int subsize = num_z_;
  if (std::abs(zbnds[1] - zbnds[0] - 1) < 1e-8) {
    full_range = true;
    --subsize;
  }
  INFO("full range " << full_range << " subsize " << subsize << " dzpn " << dzpn);
  for (int ior = 0; ior < rotator.num_orientations(); ++ior) {
    const std::vector<int>& idx = rotator.indices_[ior];
    DEBUG("stheta " << rotator.stheta_[ior]);
    scaled_relative_orientation(rotator.stheta_[ior], &s1);
    DEBUG("s1 " << s1);
    const double rh = contact.linear_interpolation(s1);
    //INFO("rh " << rh);
    const double rc = cutoff.linear_interpolation(s1);
    //INFO("rc " << rc);
    for (int idis = 0; idis < subsize; ++idis) {
      const double z = zbnds[0] + dzpn*static_cast<double>(idis);
      const double dist = rh + z*(rc - rh);
      double ener = rotator.energy(ior, dist, system);
      if (ener > NEAR_INFINITY_FLOAT) {
        ener = NEAR_INFINITY_FLOAT;
      }
      en.set_data(idx[0], idis, ener);
      if (!verbf.empty()) {
        pos.set_from_spherical({dist, rotator.stheta_[ior]});
        file << dist << "," << rotator.stheta_[ior]
          << "," << pos.str() << "," << ener << std::endl;
      }
    }
    if (full_range) {
      // Set energy to zero at cutoff (z=1)
      en.set_data(idx[0], num_z_ - 1, 0.);
      if (!verbf.empty()) {
        const int ior = subsize;
        scaled_relative_orientation(rotator.stheta_[ior], &s1);
        const double dist = cutoff.linear_interpolation(s1);
        pos.set_from_spherical({dist, rotator.stheta_[ior]});
        file << dist << "," << rotator.stheta_[subsize]
          << "," << pos.str() << "," << 0. << std::endl;
      }
    }
  }
  return en;
}

RecursiveTable5D BuildRecursiveTable::build_contact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system, const bool cutoff) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "", "");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
  RecursiveTable5D contact(argtype({{"num0", str(rotator.sizes_[0])},
                                    {"num1", str(rotator.sizes_[1])},
                                    {"num2", str(rotator.sizes_[2])},
                                    {"num3", str(rotator.sizes_[3])},
                                    {"num4", str(rotator.sizes_[4])}}));
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf);
  }
  int num_threads = -1;
//  #pragma omp parallel
  {
    num_threads = ThreadOMP().num();
  }
  std::vector<Rotator> rots(num_threads - 1);
  std::vector<System> syss(num_threads - 1);
  for (Rotator& rot : rots) rot = rotator;
  for (System& sys : syss) sys = deep_copy(*system);
//  #pragma omp parallel shared(contact)
  {
//    #pragma omp critical
    {
    Position pos(argtype({{"dimension", "3"}}));
    const int ith = ThreadOMP().thread();
    Rotator * rot;
    System * sys;
    if (ith == 0) {
      rot = &rotator;
      sys = system;
    } else {
      INFO("ith " << ith << " size " << rots.size());
      rot = &rots[ith - 1];
      sys = &syss[ith - 1];
    }
    ASSERT(rot, "er");
    ASSERT(sys, "er");
    INFO("num " << rot->num_orientations());
    for (int ior = ith; ior < rot->num_orientations(); ior += num_threads) {
      double dist;
      if (cutoff) {
        dist = rotator.cutoff_distance(ior, sys);
      } else{
        dist = rotator.contact_distance(ior, sys);
      }
      const std::vector<int>& idx = rot->indices_[ior];
      contact.set_data(idx[0], idx[1], idx[2], idx[3], idx[4], dist);
      //INFO(ior << " " << rot->contact_distance(ior, sys));
      if (!verbf.empty()) {
        //const Position& pos = system->configuration().particle(1).site(0).position();
        pos.set_from_spherical({dist, rot->stheta_[ior], rot->sphi_[ior]});
        file << dist << "," << rot->stheta_[ior] << "," << rot->sphi_[ior]
          << "," << rot->eulers_[ior].phi() << "," << rot->eulers_[ior].theta()
          << "," << rot->eulers_[ior].psi() << "," << pos.str() << std::endl;
      }
    }
    }
  }
  return contact;
}

void BuildRecursiveTable::analyze_contact_(
    const std::vector<std::vector<double> >& data, const RecursiveTable5D& contact,
    const RecursiveTable2D& contact2d, const RecursiveTable1D& contact1d,
    const std::vector<std::vector<double> >& abnd, const int dimen,
    const std::string& verbf, const std::vector<std::vector<int> >& rbins,
    double * criteria, std::vector<int> * max_bins,
    std::vector<std::vector<double> > * new_bounds,
    std::vector<double> * max_s, const bool cutoff) const {
  std::ofstream file, file2;
  if (!verbf.empty() && extra_verbose_) {
    file.open(verbf+"an");
    file2.open(verbf+"vis");
    if (dimen == 2) {
      if (static_cast<int>(data[0].size()) == 3) {
        file << "x,y" << std::endl;
      } else if (static_cast<int>(data[0].size()) == 4) {
        file << "x,y,phi" << std::endl;
      }
    } else {
      file << "x,y,z" << std::endl;
    }
  }
  std::vector<double> ms;
  std::vector<int> mb;
  Position pos, pos2;
  max_bins->clear();
  new_bounds->clear();
  *criteria = NEAR_ZERO;
  Accumulator crit_acc;
  int max_is_small = 2; // should be zero or 1
  double rh_table_max = 0.;
//  fvec5 sum;
//  fvec2 sum2d;
//  std::vector<float> sum1d;
  const bool is_2d = dimen == 2;
  const bool is_iso = is_iso_(is_2d, abnd.size());
  const Table * tab = NULL;
  if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
//    resize(contact2d.num0(), contact2d.num1(), &sum2d);
    mb.resize(2);
    tab = &contact2d;
  } else if (!is_iso && !is_2d) {
//    resize(contact.num0(), contact.num1(), contact.num2(), contact.num3(), contact.num4(), &sum);
    tab = &contact;
    mb.resize(6);
  } else if (is_iso && is_2d) {
//    sum1d.resize(contact1d.num());
    mb.resize(1);
    tab = &contact1d;
  }
  DEBUG(data.size());
  for (int idat = 0; idat < static_cast<int>(data.size()); ++idat) {
    const std::vector<double>& dat = data[idat];
    const double dist = dat[0];
    double en;
    double rh_table = -1;
    if ((is_iso && !is_2d) || (!is_iso && is_2d) ) {
      en = dat[3];
      rh_table = contact2d.linear_interpolation(dat[1], dat[2]);
    } else if (!is_iso && !is_2d) {
      en = dat[6];
      rh_table = contact.linear_interpolation(dat[1], dat[2], dat[3], dat[4], dat[5]);
    } else {
      en = dat[2];
      DEBUG("en " << en << " dat1 " << dat[1]);
      rh_table = contact1d.linear_interpolation(dat[1]);
    }

    bool overlap_predicted = false;
    if (dist < rh_table) {
      overlap_predicted = true;
    }

    bool too_small, too_big;
    if (cutoff) {
      too_small = en != 0. && !overlap_predicted;
      too_big   = en == 0. &&  overlap_predicted;
    } else {
      too_small = en > hard_limit_u_ && !overlap_predicted;
      too_big   = en < hard_limit_u_ &&  overlap_predicted;
    }
    if (!is_2d && (!verbf.empty() && extra_verbose_)) {
      pos.set_from_spherical({rh_table, (abnd[0][1]-abnd[0][0])*dat[1], (abnd[1][1]-abnd[1][0])*dat[2]});
      file << "0 " << pos.coord(0) << " " << pos.coord(1) << " " << pos.coord(2) << std::endl;
    } else if (is_2d && (!verbf.empty() && extra_verbose_)) {
      const double dtht = abnd[0][1]-abnd[0][0];
      if (static_cast<int>(dat.size()) == 3) {
        pos.set_from_spherical({rh_table, dtht*dat[1]});
        file << pos.str() << std::endl;
        pos2.set_from_spherical({dist, dtht*dat[1]});
        file2 << pos2.str() << "," << en << "," << too_small << "," << too_big << std::endl;
      } else if (static_cast<int>(dat.size()) == 4) {
        pos.set_from_spherical({rh_table, dtht*dat[1]});
        file << pos.str() << "," << (abnd[1][1] - abnd[1][0])*(dat[2] - 0.5) << std::endl;
      }
    }
    if (too_small || too_big) {
//      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
//        sum2d[contact2d.value_to_lowest_bin(0, dat[1])]
//             [contact2d.value_to_lowest_bin(1, dat[2])] = dist - rh_table;
//      } else if (is_iso && is_2d) {
//        sum1d[contact1d.value_to_lowest_bin(0, dat[1])] = dist - rh_table;
//      } else {
//        FATAL("implement");
//      }
      const double crit = std::abs(dist - rh_table);
      crit_acc.accumulate(crit);
      if (crit > *criteria) {
        DEBUG("dist " << dist << " rh_table " << rh_table << " overlap_predicted " << overlap_predicted);
        if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
          ms = {dat[1], dat[2]};
        } else if (!is_iso && !is_2d) {
          ms = {dat[1], dat[2], dat[3], dat[4], dat[5]};
        } else if (is_iso && is_2d) {
          ms = {dat[1]};
        }
        for (int rdim = 0; rdim < static_cast<int>(ms.size()); ++rdim) {
          mb[rdim] = tab->value_to_lowest_bin(rdim, ms[rdim]);
        }
        if (!find_in_list(mb, rbins)) {
          *max_s = ms;
          *criteria = crit;
          if (too_small) {
            max_is_small = true;
          } else {
            max_is_small = false;
          }
          rh_table_max = rh_table;
        }
      }
    }
  }
  INFO("max_s:" << feasst_str(*max_s) << " max_criteria:" << *criteria << " max_is_small? " << max_is_small);
  INFO("crit " << crit_acc.str());
  ASSERT(*criteria > 0, "Error");
  for (int rdim = 0; rdim < static_cast<int>(max_s->size()); ++rdim) {
    max_bins->push_back(tab->value_to_lowest_bin(rdim, (*max_s)[rdim]));
    const std::vector<double>& bn = abnd[rdim];
    const double lr = bn[0] + (bn[1] - bn[0])*tab->bin_to_value(rdim, (*max_bins)[rdim]);
    const double ur = bn[0] + (bn[1] - bn[0])*tab->bin_to_value(rdim, (*max_bins)[rdim] + 1);
    new_bounds->push_back({lr, ur});
  }
  INFO("max_bins:" << feasst_str(*max_bins));
  INFO("new_bounds:" << feasst_str(*new_bounds));
  INFO("rh_table_max:" << rh_table_max);
  if (!verbf.empty() && extra_verbose_) {
    file2 << "# " << feasst_str(*max_s) << std::endl;
  }
}

void BuildRecursiveTable::analyze_energy_(
    const std::vector<std::vector<double> >& data,
    const RecursiveTable6D& energy, const RecursiveTable3D& energy3d, const RecursiveTable2D& energy2d,
    const RecursiveTable5D& contact, const RecursiveTable2D& contact2d, const RecursiveTable1D& contact1d,
    const RecursiveTable5D& cutoff, const RecursiveTable2D& cutoff2d, const RecursiveTable1D& cutoff1d,
    const std::vector<std::vector<double> >& abnd,
    const std::string& verbf, const std::vector<std::vector<int> >& rbins,
    const bool is_2d, double * criteria,
    std::vector<int> * max_bins,
    std::vector<std::vector<double> > * new_bounds,
    std::vector<double> * new_zbnds,
    std::vector<double> * max_s) const {
  std::ofstream file;
  if (!verbf.empty() && extra_verbose_) {
    file.open(verbf+"an");
  }
  std::vector<double> ms;
  std::vector<int> mb;
  Position pos;
  max_bins->clear();
  new_bounds->clear();
  new_zbnds->clear();
  *criteria = NEAR_ZERO;
  Accumulator crit_acc;
  int max_is_small = 2; // should be zero or 1
  INFO("is 2d " << is_2d);
  const bool is_iso = is_iso_(is_2d, abnd.size());
  INFO("is iso " << is_iso);
  const Table * tab;
  if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
    tab = &energy3d;
    mb.resize(3);
  } else if (!is_iso && !is_2d) {
    tab = &energy;
    mb.resize(6);
  } else if (is_iso && is_2d) {
    tab = &energy2d;
    mb.resize(2);
  } else { // if (!is_iso && is_2d) {
    FATAL("implement");
  }
  INFO(data.size());
  for (int idat = 0; idat < static_cast<int>(data.size()); ++idat) {
    const std::vector<double>& dat = data[idat];
    const double dist = dat[0];
    double rh, rc;
    if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
      rh = contact2d.linear_interpolation(dat[1], dat[2]);
      rc = cutoff2d.linear_interpolation(dat[1], dat[2]);
    } else if (!is_iso && !is_2d) {
      rh = contact.linear_interpolation(dat[1], dat[2], dat[3], dat[4], dat[5]);
      rc = cutoff.linear_interpolation(dat[1], dat[2], dat[3], dat[4], dat[5]);
    } else if (is_iso && is_2d) {
      rh = contact1d.linear_interpolation(dat[1]);
      rc = cutoff1d.linear_interpolation(dat[1]);
    } else {
      FATAL("implement");
    }
    const double z = (dist - rh)/(rc - rh);
    //INFO("dist " << dist << " rh " << rh << " rc " << rc);
    if (z >= 0. && z <= 1.) {
      double en;
      double en_table = -1;
      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
        en = dat[3];
        en_table = energy3d.linear_interpolation(dat[1], dat[2], z);
      } else if (!is_iso && !is_2d) {
        en = dat[6];
        en_table = energy.linear_interpolation(dat[1], dat[2], dat[3], dat[4], dat[5], z);
      } else if (is_iso && is_2d) {
        en = dat[2];
        en_table = energy2d.linear_interpolation(dat[1], z);
      } else {
        FATAL("implement");
      }
      if (!verbf.empty() && extra_verbose_) {
        if (is_2d) {
          pos.set_from_spherical({en_table, PI*(2*dat[1]-1)});
          file << pos.str() << std::endl;
        } else {
          pos.set_from_spherical({en_table, PI*(2*dat[1]-1), PI*dat[2]});
          file << "0 " << pos.coord(0) << " " << pos.coord(1) << " " << pos.coord(2) << std::endl;
          //file << pos.str() << std::endl;
        }
      }
      const double crit = std::abs(std::exp(-beta_*en) - std::exp(-beta_*en_table));
      //const double crit = dist*dist*std::abs(std::exp(-beta_*en) - std::exp(-beta_*en_table));
      crit_acc.accumulate(crit);
      if (crit > *criteria) {
        if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
          ms = {dat[1], dat[2], z};
        } else if (!is_iso && !is_2d) {
          ms = {dat[1], dat[2], dat[3], dat[4], dat[5], z};
        } else if (is_iso && is_2d) {
          ms = {dat[1], z};
        } else if (!is_iso && is_2d) {
          FATAL("implement");
        }
        for (int rdim = 0; rdim < static_cast<int>(ms.size()); ++rdim) {
          mb[rdim] = tab->value_to_lowest_bin(rdim, ms[rdim]);
        }
        if (!find_in_list(mb, rbins)) {
          *max_s = ms;
          *criteria = crit;
        }
        //INFO("en " << en << " en_table " << en_table << " dist " << dist);
      }
    }
  }
  INFO("max_s:" << feasst_str(*max_s) << " max_criteria:" << *criteria << " max_is_small? " << max_is_small);
  INFO("crit " << crit_acc.str());
  ASSERT(*criteria > 0, "Error");
  int rdim;
  for (rdim = 0; rdim < static_cast<int>(max_s->size() - 1); ++rdim) {
    max_bins->push_back(tab->value_to_lowest_bin(rdim, (*max_s)[rdim]));
    const std::vector<double>& bn = abnd[rdim];
    const double lr = bn[0] + (bn[1] - bn[0])*tab->bin_to_value(rdim, (*max_bins)[rdim]);
    const double ur = bn[0] + (bn[1] - bn[0])*tab->bin_to_value(rdim, (*max_bins)[rdim] + 1);
    new_bounds->push_back({lr, ur});
  }
  //INFO("rdim " << rdim);
  max_bins->push_back(tab->value_to_lowest_bin(rdim, (*max_s)[rdim]));
  const double lz = tab->bin_to_value(rdim, (*max_bins)[rdim]);
  const double uz = tab->bin_to_value(rdim, (*max_bins)[rdim] + 1);
  *new_zbnds = {lz, uz};
  INFO("max_bins:" << feasst_str(*max_bins));
  INFO("new_bounds:" << feasst_str(*new_bounds));
  INFO("new_zbnds:" << feasst_str(*new_zbnds));
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

void BuildRecursiveTable::run(MonteCarlo * mc) {
  std::vector<std::vector<double> > data;
  read_data_(mayer_training_file_, &data);
//  std::vector<bool> datakeep(data.size(), true);
  const double lower = global_lower_(hard_limit_u_, &data);
  const double upper = global_upper_(cutoff_, &data);
  std::stringstream ss;
  // Todo for multi-site tableing (5/6D , not 1D data)
  // Separate contact table necessary? Separate cutoff table? (in trimer case, cutoff varies, in all atom protein, maybe not)
  // Is a separate MayerSampling simulation for hard particle interactions to determine rh necessary, or can you simply use the full potential for training. Contact only might be an easier place to test? But will its varied resolution mesh with a separate table? z will change with rh, so rh/rc has to be determined first and then not changed. Set hard spheres at the u=hardlimit distance, define rh, test against hard sphere b2 for rh-only excluded volume
  // Set bounds for all dimensions, used with build_table and analyze_table
  // Use Rotator, but might need to accept custom bounds
  // Parallelize
  double criteria;
  int iteration = 0;
  std::string verbf = verbose_name_(iteration);
  std::vector<std::vector<double> > new_bounds;
  std::vector<int> max_bins;
  std::vector<double> max_s;
  if (data[0].size() == 2) {
    RecursiveTable1D energy_table(argtype({{"num", str(num_z_)}}));
    build_table_({{lower, upper}}, data, &energy_table, mc);
    int max_bin = analyze_table_(lower, upper, data, energy_table, beta_, verbf, &criteria);
    double last_criteria = criteria;
    RecursiveTable1D nested(argtype({{"num", str(num_z_)}}));
    std::vector<int> bins = {max_bin};
    const double zfac = upper - lower;
    while (criteria > min_criteria_) {
      ++iteration;
      verbf = verbose_name_(iteration);
      const double lr = lower + zfac*energy_table.bin_value(max_bin);
      const double ur = lower + zfac*energy_table.bin_value(max_bin+1);
      build_table_({{lr, ur}}, data, &nested, mc);
      energy_table.insert(max_bin, nested);
      max_bin = analyze_table_(lower, upper, data, energy_table, beta_, verbf, &criteria);
      ASSERT(criteria < last_criteria, "Adding resolution made crit:" << criteria << " worse than last:" << last_criteria);
      last_criteria = criteria;
      ASSERT(iteration < 1e4, "Cannot get criteria:" << criteria <<
        " below min_criteria:" << min_criteria_ << " within iteration:"
        << iteration);
      if (criteria > min_criteria_) {
        if (find_in_list(max_bin, bins)) {
          INFO("At iteration:" << iteration <<
            ", max_bin:" << max_bin << " was already nested but found to have "
            << "criteria:" << criteria << ", which is higher than min_criteria:"
            << min_criteria_);
          break;
        }
      }
      bins.push_back(max_bin);
    }
    ModelRecursiveTable pot;
    pot.lower_ = {{lower}};
    pot.upper_ = {{upper}};
    pot.energy_ = {{energy_table}};
    pot.serialize(ss);
  } else if (data[0].size() == 3 || data[0].size() == 4 || data[0].size() == 7
          || data[0].size() == 10 || data[0].size() == 12) {
    // for anisotropic particles, contact and cutoff must be build first, then energy which is based on z=(r-rh)/(rc-rh)
    // the resolution does not have to be the same for any of these tables
    const Configuration& config = mc->system().configuration();
    const std::vector<std::vector<double> >& abnd = Rotator().gen_global_bounds(config);
    DEBUG("abnd " << feasst_str(abnd) << " " << abnd.size());
    const bool is_2d  = static_cast<int>(config.dimension()) == 2;
    DEBUG("is_2d: " << is_2d);
    const bool is_iso = is_iso_(is_2d, abnd.size());
    RecursiveTable1D contact1d;
    RecursiveTable2D contact2d;
    RecursiveTable5D contact;
    if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
      ASSERT(data[0].size() == 4, data[0].size());
      contact2d = build_2dcontact_(abnd, verbf, mc->get_system());
      INFO("num_data: " << contact2d.num_data());
      INFO("percent_nested: " << contact2d.percent_nested());
    } else if (is_iso && is_2d) {
      ASSERT(data[0].size() == 3, data[0].size());
      contact1d = build_1dcontact_(abnd, verbf, mc->get_system());
      INFO("percent_nested: " << contact1d.percent_nested());
    } else if (!is_iso && !is_2d) {
      contact = build_contact_(abnd, verbf, mc->get_system());
    }
    std::vector<std::vector<int> > rbins;
    analyze_contact_(data, contact, contact2d, contact1d, abnd, config.dimension(), verbf, rbins, &criteria, &max_bins, &new_bounds, &max_s);
    double last_criteria = criteria;
    INFO("max_s " << feasst_str(max_s));
    rbins = {max_bins};
    INFO("criteria: " << criteria);
    //while (false) {
    RecursiveTable5D nested;
    RecursiveTable2D nested2d;
    RecursiveTable1D nested1d;
    //WARN("fixed 1 iteration");
    //while (iteration <= 1 && min_criteria_ < 200) {
    //while (iteration == 0 && min_criteria_ < 200) {
    //while (false) {
    while (criteria > min_criteria_) {
      ++iteration;
      verbf = verbose_name_(iteration);
      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
        nested2d = build_2dcontact_(new_bounds, verbf, mc->get_system());
        INFO("adding recursive table at " << feasst_str(max_bins));
        contact2d.insert(max_bins[0], max_bins[1], nested2d);
        INFO("num_data: " << contact2d.num_data());
        INFO("percent_nested: " << contact2d.percent_nested());
      } else if (!is_iso && !is_2d) {
        nested = build_contact_(new_bounds, verbf, mc->get_system());
        INFO("adding recursive table at " << feasst_str(max_bins));
        contact.insert(max_bins[0], max_bins[1], max_bins[2], max_bins[3], max_bins[4], nested);
        INFO("percent_nested: " << contact.percent_nested());
      } else if (is_iso && is_2d) {
        nested1d = build_1dcontact_(new_bounds, verbf, mc->get_system());
        INFO("adding recursive table at " << feasst_str(max_bins));
        contact1d.insert(max_bins[0], nested1d);
        INFO("percent_nested: " << contact1d.percent_nested());
      }
      analyze_contact_(data, contact, contact2d, contact1d, abnd, config.dimension(), verbf, rbins, &criteria, &max_bins, &new_bounds, &max_s);
      if (!is_2d) {
        WARN("no need for max_s, just testing");
        Position pos(argtype({{"dimension", "3"}}));
        pos.set_from_spherical({1, PI*(2*max_s[0]-1), PI*max_s[1]});
        INFO("max_s " << feasst_str(max_s) << " cart: " << pos.str());
        if (is_iso) {
          const double rh_table = contact2d.linear_interpolation(max_s[0], max_s[1]);
          INFO("testing new rh_table:" << rh_table);
        } else {
          const double rh_table = contact.linear_interpolation(max_s[0], max_s[1], max_s[2], max_s[3], max_s[4]);
          INFO("testing new rh_table:" << rh_table);
        }
      }
      ASSERT(criteria <= last_criteria, "Adding resolution made crit:" << criteria << " worse than last:" << last_criteria);
      last_criteria = criteria;
      ASSERT(iteration < 1e4, "Cannot get criteria:" << criteria <<
        " below min_criteria:" << min_criteria_ << " within iteration:"
        << iteration);
      if (criteria > min_criteria_) {
        if (find_in_list(max_bins, rbins)) {
          INFO("At iteration:" << iteration <<
            ", max_bins:" << feasst_str(max_bins) << " was already nested but found to have "
            << "criteria:" << criteria << ", which is higher than min_criteria:"
            << min_criteria_);
          break;
        }
      }
      rbins.push_back(max_bins);
      INFO("criteria: " << criteria);
    }

    // write to file
    RecursiveTable vis;
    INFO("vis ignore en? " << vis.ignore_energy());
    if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
      vis.contact2d_ = {{contact2d}};
    } else if (!is_iso && !is_2d) {
      vis.contact_ = {{contact}};
    } else if (is_iso && is_2d) {
      vis.contact1d_ = {{contact1d}};
    }

    if (!contact_only_) {
      iteration = 0;
      verbf = verbose_name_(iteration) + "cut";
      INFO("building cutoff");
      RecursiveTable1D cutoff1d;
      RecursiveTable2D cutoff2d;
      RecursiveTable5D cutoff;
      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
        cutoff2d = build_2dcontact_(abnd, verbf, mc->get_system(), true);
      } else if (!is_iso && !is_2d) {
        cutoff = build_contact_(abnd, verbf, mc->get_system(), true);
      } else if (is_iso && is_2d) {
        cutoff1d = build_1dcontact_(abnd, verbf, mc->get_system(), true);
      }
      rbins.clear();
      analyze_contact_(data, cutoff, cutoff2d, cutoff1d, abnd, config.dimension(), verbf, rbins, &criteria, &max_bins, &new_bounds, &max_s, true);
      last_criteria = criteria;
      rbins = {max_bins};
      //while (false) {
      while (criteria > min_criteria_) {
        ++iteration;
        verbf = verbose_name_(iteration) + "cut";
        if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
          nested2d = build_2dcontact_(new_bounds, verbf, mc->get_system(), true);
          cutoff2d.insert(max_bins[0], max_bins[1], nested2d);
          INFO("percent_nested: " << cutoff2d.percent_nested());
        } else if (!is_iso && !is_2d) {
          nested = build_contact_(new_bounds, verbf, mc->get_system(), true);
          cutoff.insert(max_bins[0], max_bins[1], max_bins[2], max_bins[3], max_bins[4], nested);
          INFO("percent_nested: " << cutoff.percent_nested());
        } else if (is_iso && is_2d) {
          nested1d = build_1dcontact_(new_bounds, verbf, mc->get_system(), true);
          cutoff1d.insert(max_bins[0], nested1d);
          INFO("percent_nested: " << cutoff1d.percent_nested());
        } else {
          FATAL("implement");
        }
        analyze_contact_(data, cutoff, cutoff2d, cutoff1d, abnd, config.dimension(), verbf, rbins, &criteria, &max_bins, &new_bounds, &max_s, true);
        ASSERT(criteria <= last_criteria, "Adding resolution made crit:" << criteria << " worse than last:" << last_criteria);
        last_criteria = criteria;
        ASSERT(iteration < 1e4, "Cannot get criteria:" << criteria <<
          " below min_criteria:" << min_criteria_ << " within iteration:"
          << iteration);
        if (criteria > min_criteria_) {
          if (find_in_list(max_bins, rbins)) {
            INFO("At iteration:" << iteration <<
              ", max_bins:" << feasst_str(max_bins) << " was already nested but found to have "
              << "criteria:" << criteria << ", which is higher than min_criteria:"
              << min_criteria_);
            break;
          }
        }
        rbins.push_back(max_bins);
        INFO("criteria: " << criteria);
      }
      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
        vis.cutoff2d_ = {{cutoff2d}};
      } else if (!is_iso && !is_2d) {
        vis.cutoff_ = {{cutoff}};
      } else if (is_iso && is_2d) {
        vis.cutoff1d_ = {{cutoff1d}};
      } else {
        FATAL("implement");
      }
      INFO("building energy");
      iteration = 0;
      verbf = verbose_name_(iteration) + "en";
      RecursiveTable2D energy2d;
      RecursiveTable3D energy3d;
      RecursiveTable6D energy;
      RecursiveTable2D en_nested2d;
      RecursiveTable3D en_nested3d;
      RecursiveTable6D en_nested;
      std::vector<double> new_zbnds;
      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
        energy3d = build_3denergy_(abnd, {0., 1.}, verbf, contact2d, cutoff2d, mc->get_system());
      } else if (!is_iso && !is_2d) {
        energy = build_energy_(abnd, {0., 1.}, verbf, contact, cutoff, mc->get_system());
      } else if (is_iso && is_2d) {
        energy2d = build_2denergy_(abnd, {0., 1.}, verbf, contact1d, cutoff1d, mc->get_system());
      } else {
        FATAL("implement");
      }
      rbins.clear();
      analyze_energy_(data, energy, energy3d, energy2d, contact, contact2d, contact1d, cutoff, cutoff2d, cutoff1d, abnd, verbf, rbins, is_2d, &criteria, &max_bins, &new_bounds, &new_zbnds, &max_s);
      last_criteria = criteria;
      rbins = {max_bins};
      while (criteria > min_criteria_energy_) {
        ++iteration;
        verbf = verbose_name_(iteration) + "en";
        INFO("is_iso " << is_iso);
        //if (is_iso && !is_2d) {
        if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
          en_nested3d = build_3denergy_(new_bounds, new_zbnds, verbf, contact2d, cutoff2d, mc->get_system());
          //INFO("mb " << feasst_str(max_bins));
          energy3d.insert(max_bins[0], max_bins[1], max_bins[2], en_nested3d);
          INFO("percent_nested: " << energy3d.percent_nested());
        } else if (!is_iso && !is_2d) {
          en_nested = build_energy_(new_bounds, new_zbnds, verbf, contact, cutoff, mc->get_system());
          energy.insert(max_bins[0], max_bins[1], max_bins[2], max_bins[3], max_bins[4], max_bins[5], en_nested);
          INFO("percent_nested: " << energy.percent_nested());
        } else if (is_iso && is_2d) {
          en_nested2d = build_2denergy_(new_bounds, new_zbnds, verbf, contact1d, cutoff1d, mc->get_system());
          energy2d.insert(max_bins[0], max_bins[1], en_nested2d);
          INFO("percent_nested: " << energy2d.percent_nested());
        } else {
          FATAL("implement");
        }
        analyze_energy_(data, energy, energy3d, energy2d, contact, contact2d, contact1d, cutoff, cutoff2d, cutoff1d, abnd, verbf, rbins, is_2d, &criteria, &max_bins, &new_bounds, &new_zbnds, &max_s);
        INFO("criteria " << criteria << " last_criteria " << last_criteria);
        ASSERT(criteria <= last_criteria, "Adding resolution made crit:" << criteria << " worse than last:" << last_criteria);
        last_criteria = criteria;
        ASSERT(iteration < 1e4, "Cannot get criteria:" << criteria <<
          " below min_criteria:" << min_criteria_energy_ << " within iteration:"
          << iteration);
        if (criteria > min_criteria_energy_) {
          if (find_in_list(max_bins, rbins)) {
            INFO("At iteration:" << iteration <<
              ", max_bins:" << feasst_str(max_bins) << " was already nested but found to have "
              << "criteria:" << criteria << ", which is higher than min_criteria:"
              << min_criteria_energy_);
            break;
          }
        }
        rbins.push_back(max_bins);
        INFO("criteria: " << criteria);
      }
      if ( (is_iso && !is_2d) || (!is_iso && is_2d) ) {
        vis.energy3d_ = {{energy3d}};
      } else if (!is_iso && !is_2d) {
        vis.energy_ = {{energy}};
      } else if (is_iso && is_2d) {
        vis.energy2d_ = {{energy2d}};
      } else {
        FATAL("implement");
      }
    }
    INFO("vis ignore en? " << vis.ignore_energy());
    vis.serialize(ss);
  } else {
    FATAL("unrecognized data size:" << data[0].size());
  }

  // write serialization to file
  ASSERT(!output_file_.empty(), "Error");
  std::ofstream file2(output_file_);
  ASSERT(file2.good(), "Error");
  file2 << ss.str();
  data.clear();
}

}  // namespace feasst
