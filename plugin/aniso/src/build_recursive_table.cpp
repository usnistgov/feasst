#include <fstream>
#include "threads/include/thread_omp.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
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
#include "aniso/include/visit_model_inner_recursive_table.h"
#include "aniso/include/build_recursive_table.h"

namespace feasst {

BuildRecursiveTable::BuildRecursiveTable(argtype * args) {
  class_name_ = "BuildRecursiveTable";
  mayer_training_file_ = str("mayer_training_file", args, "");
  cutoff_ = dble("cutoff", args, -1);
  output_file_ = str("output_file", args, "");
  verbose_file_ = str("verbose_file", args, "");
  hard_limit_u_ = dble("hard_limit_u", args, 100);
  size_ = dble("size", args, 5);
  num_orientations_per_pi_ = dble("num_orientations_per_pi", args, 5);
  beta_ = dble("beta", args, 1.);
  min_criteria_ = dble("min_criteria", args, 0.03);
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
  feasst_deserialize(&size_, istr);
  feasst_deserialize(&num_orientations_per_pi_, istr);
  feasst_deserialize(&beta_, istr);
  feasst_deserialize(&min_criteria_, istr);
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
  feasst_serialize(size_, ostr);
  feasst_serialize(num_orientations_per_pi_, ostr);
  feasst_serialize(beta_, ostr);
  feasst_serialize(min_criteria_, ostr);
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
  const int max_bin = table.value_to_lowest_bin(max_z);
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

RecursiveTable2D BuildRecursiveTable::build_2dcontact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "tr6d_fixed", "tr6d_contact");
  rotator.gen_unique_orientations(num_orientations_per_pi_, system, bounds);
  RecursiveTable2D contact(argtype({{"num0", str(rotator.sizes_[0])},
                                    {"num1", str(rotator.sizes_[1])}}));
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
      const double dist = rot->contact_distance(ior, sys);
      const std::vector<int>& idx = rot->indices_[ior];
      contact.set_data(idx[0], idx[1], dist);
      //INFO(ior << " " << rot->contact_distance(ior, sys));
      if (!verbf.empty()) {
        //const Position& pos = system->configuration().particle(1).site(0).position();
        pos.set_from_spherical({dist, rot->stheta_[ior], rot->sphi_[ior]});
        file << dist << "," << rot->stheta_[ior] << "," << rot->sphi_[ior]
          << "," << pos.str() << std::endl;
      }
    }
    }
  }
  return contact;
}

RecursiveTable5D BuildRecursiveTable::build_contact_(const std::vector<std::vector<double> >& bounds, const std::string& verbf, System * system) {
  Rotator rotator({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u_)}});
  rotator.init(system, "tr6d_fixed", "tr6d_contact");
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
      const double dist = rot->contact_distance(ior, sys);
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
    const std::vector<std::vector<double> >& data, const RecursiveTable5D& contact, const RecursiveTable2D& contact2d,
    const std::vector<std::vector<double> >& abnd, const std::string& verbf, double * criteria,
    std::vector<int> * max_bins,
    std::vector<std::vector<double> > * new_bounds,
    std::vector<double> * max_s) const {
  std::ofstream file;
  if (!verbf.empty()) {
    file.open(verbf+"an");
  }
  Position pos;
  max_bins->clear();
  new_bounds->clear();
  *criteria = -1.;
  int max_is_small = 2; // should be zero or 1
  double rh_table_max = 0.;
  fvec5 sum;
  fvec2 sum2d;
  const bool is_2d = static_cast<int>(abnd.size()) == 2;
  if (is_2d) {
    resize(contact2d.num0(), contact2d.num1(), &sum2d);
  } else {
    resize(contact.num0(), contact.num1(), contact.num2(), contact.num3(), contact.num4(), &sum);
  }
  for (int idat = 0; idat < static_cast<int>(data.size()); ++idat) {
    const std::vector<double>& dat = data[idat];
    const double dist = dat[0];
    double en;
    double rh_table = -1;
    if (is_2d) {
      en = dat[3];
      rh_table = contact2d.linear_interpolation(dat[1], dat[2]);
    } else {
      en = dat[6];
      rh_table = contact.linear_interpolation(dat[1], dat[2], dat[3], dat[4], dat[5]);
    }
    if (!verbf.empty()) {
      // WARN: below assumes i != j symmetry
      pos.set_from_spherical({rh_table, PI*(2*dat[1]-1), PI*dat[2]});
      file << "0 " << pos.coord(0) << " " << pos.coord(1) << " " << pos.coord(2) << std::endl;
      //file << pos.str() << std::endl;
    }
    bool overlap_predicted = false;
    if (dist < rh_table) {
      overlap_predicted = true;
    }
    const bool too_small = en > 1 && !overlap_predicted;
    const bool too_big   = en < 1 &&  overlap_predicted;
    if (too_small || too_big) {
      if (is_2d) {
        sum2d[contact2d.value_to_lowest_bin(0, dat[1])]
             [contact2d.value_to_lowest_bin(1, dat[2])] = dist - rh_table;
      } else {
        FATAL("implement");
      }
      const double crit = std::abs(dist - rh_table);
      if (crit > *criteria) {
        if (is_2d) {
          *max_s = {dat[1], dat[2]};
        } else {
          *max_s = {dat[1], dat[2], dat[3], dat[4], dat[5]};
        }
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
  ASSERT(*criteria > 0, "Error");
  INFO("max_s:" << feasst_str(*max_s) << " max_criteria:" << *criteria << " max_is_small? " << max_is_small);
  if (is_2d) {
//    //float max
//    *criteria = std::numeric_limits<float>::min();
//    for (int i0 = 0; i0 < static_cast<int>(sum2d.size()); ++i0) {
//      for (int i1 = 0; i1 < static_cast<int>(sum2d[i0].size()); ++i1) {
//        const double abssum2d = std::abs(sum2d[i0][i1]);
//        if (*criteria < abssum2d) {
//          *criteria = abssum2d;
//          *max_s = {contact2d.bin_to_value(0, i0), contact2d.bin_to_value(1, i1)};
//        }
//      }
//    }

    for (int rdim = 0; rdim < static_cast<int>(max_s->size()); ++rdim) {
      max_bins->push_back(contact2d.value_to_lowest_bin(rdim, (*max_s)[rdim]));
      const std::vector<double>& bn = abnd[rdim];
      const double lr = bn[0] + (bn[1] - bn[0])*contact2d.bin_to_value(rdim, (*max_bins)[rdim]);
      const double ur = bn[0] + (bn[1] - bn[0])*contact2d.bin_to_value(rdim, (*max_bins)[rdim] + 1);
      new_bounds->push_back({lr, ur});
    }
  } else {
    for (int rdim = 0; rdim < static_cast<int>(max_s->size()); ++rdim) {
      max_bins->push_back(contact.value_to_lowest_bin(rdim, (*max_s)[rdim]));
      const std::vector<double>& bn = abnd[rdim];
      const double lr = bn[0] + (bn[1] - bn[0])*contact.bin_to_value(rdim, (*max_bins)[rdim]);
      const double ur = bn[0] + (bn[1] - bn[0])*contact.bin_to_value(rdim, (*max_bins)[rdim] + 1);
      new_bounds->push_back({lr, ur});
    }
  }
  INFO("max_bins:" << feasst_str(*max_bins));
  INFO("new_bounds:" << feasst_str(*new_bounds));
  INFO("rh_table_max:" << rh_table_max);
}

void BuildRecursiveTable::run(MonteCarlo * mc) {
  std::vector<std::vector<double> > data;
  read_data_(mayer_training_file_, &data);
  const double lower = global_lower_(hard_limit_u_, &data);
  const double upper = global_upper_(cutoff_, &data);
  std::stringstream ss;
  // Todo for multi-site tableing (5/6D , not 1D data)
  // Separate contact table necessary? Separate cutoff table? (in trimer case, cutoff varies, in all atom protein, maybe not)
  // Is a separate MayerSampling simulation for hard particle interactions to determine rh necessary, or can you simply use the full potential for training. Contact only might be an easier place to test? But will its varied resolution mesh with a separate table? z will change with rh, so rh/rc has to be determined first and then not changed. Set hard spheres at the u=hardlimit distance, define rh, test against hard sphere b2 for rh-only excluded volume
  // Set bounds for all dimensions, used with build_table and analyze_table
  // Use Rotator, but might need to accept custom bounds
  // Parallelize
  double criteria, last_criteria = NEAR_INFINITY;
  int iteration = 0;
  std::string verbf = verbose_name_(iteration);
  std::vector<std::vector<double> > new_bounds;
  std::vector<int> max_bins;
  std::vector<double> max_s;
  if (data[0].size() == 2) {
    RecursiveTable1D energy_table(argtype({{"num", str(size_)}}));
    build_table_({{lower, upper}}, data, &energy_table, mc);
    int max_bin = analyze_table_(lower, upper, data, energy_table, beta_, verbf, &criteria);
    RecursiveTable1D nested(argtype({{"num", str(size_)}}));
    std::vector<int> bins = {max_bin};
    const double zfac = upper - lower;
    while (criteria > min_criteria_) {
      ++iteration;
      verbf = verbose_name_(iteration);
      const double lr = lower + zfac*energy_table.bin_to_value(max_bin);
      const double ur = lower + zfac*energy_table.bin_to_value(max_bin+1);
      build_table_({{lr, ur}}, data, &nested, mc);
      energy_table.insert(max_bin, nested);
      max_bin = analyze_table_(lower, upper, data, energy_table, beta_, verbf, &criteria);
      ASSERT(criteria < last_criteria, "Adding resolution made crit:" << criteria << " worse than last:" << last_criteria);
      last_criteria = criteria;
      ASSERT(iteration < 1e4, "Cannot get criteria:" << criteria <<
        " below min_criteria:" << min_criteria_ << " within iteration:"
        << iteration);
      if (criteria > min_criteria_) {
        ASSERT(!find_in_list(max_bin, bins), "At iteration:" << iteration <<
          ", max_bin:" << max_bin << " was already nested but found to have "
          << "criteria:" << criteria << ", which is higher than min_criteria:"
          << min_criteria_);
      }
      bins.push_back(max_bin);
    }
    ModelRecursiveTable pot;
    pot.lower_ = {{lower}};
    pot.upper_ = {{upper}};
    pot.energy_ = {{energy_table}};
    pot.serialize(ss);
//  } else if (data[4].size() == 4) {
//    const std::vector<std::vector<double> >& abnd = Rotator().gen_global_bounds(mc->system().configuration());
//    INFO("abnd " << feasst_str(abnd));
//    RecursiveTable5D contact;
//    analyze_contact_(data, contact, contact2d, abnd, verbf, &criteria, &max_bins, &new_bounds, &max_s);
//    std::vector<std::vector<int> > rbins = {max_bins};
  } else if (data[4].size() == 4 || data[0].size() == 7 || data[0].size() == 10 || data[0].size() == 12) {
    const std::vector<std::vector<double> >& abnd = Rotator().gen_global_bounds(mc->system().configuration());
    INFO("abnd " << feasst_str(abnd));
    bool is_2d = static_cast<int>(abnd.size()) == 2;
    RecursiveTable2D contact2d;
    RecursiveTable5D contact;
    if (is_2d) {
      contact2d = build_2dcontact_(abnd, verbf, mc->get_system());
      INFO("num_data: " << contact2d.num_data());
      INFO("percent_nested: " << contact2d.percent_nested());
    } else {
      contact = build_contact_(abnd, verbf, mc->get_system());
    }
    analyze_contact_(data, contact, contact2d, abnd, verbf, &criteria, &max_bins, &new_bounds, &max_s);
    INFO("max_s " << feasst_str(max_s));
    std::vector<std::vector<int> > rbins = {max_bins};
    INFO("criteria: " << criteria);
    //while (false) {
    RecursiveTable5D nested;
    RecursiveTable2D nested2d;
    while (criteria > min_criteria_) {
      ++iteration;
      verbf = verbose_name_(iteration);
      if (is_2d) {
        nested2d = build_2dcontact_(new_bounds, verbf, mc->get_system());
        INFO("adding recursive table at " << feasst_str(max_bins));
        contact2d.insert(max_bins[0], max_bins[1], nested2d);
        INFO("num_data: " << contact2d.num_data());
        INFO("percent_nested: " << contact2d.percent_nested());
      } else {
        nested = build_contact_(new_bounds, verbf, mc->get_system());
        INFO("adding recursive table at " << feasst_str(max_bins));
        contact.insert(max_bins[0], max_bins[1], max_bins[2], max_bins[3], max_bins[4], nested);
      }
      analyze_contact_(data, contact, contact2d, abnd, verbf, &criteria, &max_bins, &new_bounds, &max_s);
      { WARN("no need for max_s, just testing");
        Position pos(argtype({{"dimension", "3"}}));
        pos.set_from_spherical({1, PI*(2*max_s[0]-1), PI*max_s[1]});
        INFO("max_s " << feasst_str(max_s) << " cart: " << pos.str());
        if (is_2d) {
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
        ASSERT(!find_in_list(max_bins, rbins), "At iteration:" << iteration <<
          ", max_bins:" << feasst_str(max_bins) << " was already nested but found to have "
          << "criteria:" << criteria << ", which is higher than min_criteria:"
          << min_criteria_);
      }
      rbins.push_back(max_bins);
      INFO("criteria: " << criteria);
    }

    // write to file
    VisitModelInnerRecursiveTable vis;
    if (is_2d) {
      vis.contact2d_ = {{contact2d}};
    } else {
      vis.contact_ = {{contact}};
    }
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
