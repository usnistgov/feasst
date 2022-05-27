#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "system/include/energy_map.h"

namespace feasst {

std::map<std::string, std::shared_ptr<EnergyMap> >& EnergyMap::deserialize_map() {
  static std::map<std::string, std::shared_ptr<EnergyMap> >* ans =
     new std::map<std::string, std::shared_ptr<EnergyMap> >();
  return *ans;
}

EnergyMap::EnergyMap(argtype * args) {
  default_value_ = dble("default_value", args, 0.);
}
EnergyMap::EnergyMap(argtype args) : EnergyMap(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void EnergyMap::clear(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  resize_(part1_index, site1_index, part2_index, site2_index);
  *smap_new_(part1_index, site1_index, part2_index, site2_index) =
    std::vector<double>(2 + dimen_, default_value());
  *smap_new_(part2_index, site2_index, part1_index, site1_index) =
    std::vector<double>(2 + dimen_, default_value());
}

double EnergyMap::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int site1_type,
    const int part2_index,
    const int site2_index,
    const int site2_type,
    const double squared_distance,
    const Position * pbc) {
  resize_(part1_index, site1_index, part2_index, site2_index);
  std::vector<double> * smap1 =
    smap_new_(part1_index, site1_index, part2_index, site2_index);
  (*smap1)[0] = energy;
  (*smap1)[1] = squared_distance;
  if (pbc->dimension() > 0) {
    for (int dim = 0; dim < dimen_; ++dim) {
      (*smap1)[2 + dim] = pbc->coord(dim);
    }
  }
  std::vector<double> * smap2 =
    smap_new_(part2_index, site2_index, part1_index, site1_index);
  (*smap2)[0] = energy;
  (*smap2)[1] = squared_distance;
  if (pbc->dimension() > 0) {
    for (int dim = 0; dim < dimen_; ++dim) {
      (*smap2)[2 + dim] = -1.*pbc->coord(dim);
    }
  }
  return energy;
}

void EnergyMap::precompute(Configuration * config) {
  site_max_ = config->max_sites_in_any_particle();
  DEBUG("site_max_ " << site_max_);
  dimen_ = config->dimension();
}

double EnergyMap::total_energy() const {
  double en = 0;
  //DEBUG(map().size() << " " << map()[0].size() << " " << map()[0][0].size() << " " << map()[0][0][0].size() << " " << map()[0][0][0][0].size());
  for (const std::vector<std::vector<std::vector<std::vector<double> > > >& vec4 : map()) {
    for (const std::vector<std::vector<std::vector<double> > >& vec3 : vec4) {
      for (const std::vector<std::vector<double> >& vec2 : vec3) {
        for (const std::vector<double>& vec1 : vec2) {
          //DEBUG(vec1.size());
          en += vec1[0];
        }
      }
    }
  }
  return 0.5*en;
}

std::shared_ptr<EnergyMap> EnergyMap::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void EnergyMap::serialize_energy_map_(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(945, ostr);
  feasst_serialize(default_value_, ostr);
  feasst_serialize(site_max_, ostr);
  feasst_serialize(dimen_, ostr);
  feasst_serialize_fstobj(data_, ostr);
}

EnergyMap::EnergyMap(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 945, "version mismatch: " << version);
  feasst_deserialize(&default_value_, istr);
  feasst_deserialize(&site_max_, istr);
  feasst_deserialize(&dimen_, istr);
  feasst_deserialize_fstobj(&data_, istr);
}

void EnergyMap::select_cluster(const NeighborCriteria& cluster_criteria,
                              const Configuration& config,
                              const int particle_node,
                              Select * cluster,
                              const Position& frame_of_reference) const {
  FATAL("not implemented");
}

bool EnergyMap::is_cluster_changed(const NeighborCriteria& cluster_criteria,
    const Select& select,
    const Configuration& config) const {
  FATAL("not implemented");
}

const NeighborCriteria& EnergyMap::neighbor_criteria() const {
  FATAL("not implemented");
}

void EnergyMap::neighbors(
    const NeighborCriteria& neighbor_criteria,
    const Configuration& config,
    const int target_particle,
    const int target_site,
    const int given_site_index,
    Select * neighbors,
    const int new_map) const {
  FATAL("not implemented");
}

void EnergyMap::synchronize_(const EnergyMap& map, const Select& perturbed) {
  FATAL("not implemented");
}

double EnergyMap::energy(const int part1_index, const int site1_index) const {
  FATAL("not implemented");
}

//const std::vector<double>& EnergyMap::map(const int part1, const int part2,
//    const int site1, const int site2) const {
//  FATAL("not implemented");
//}

std::shared_ptr<EnergyMap> EnergyMap::create(std::istream& istr) const {
  FATAL("not implemented");
}
std::shared_ptr<EnergyMap> EnergyMap::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<EnergyMap> EnergyMap::factory(const std::string name, argtype * args) {
  DEBUG("name: " << name << ", args: " << str(*args));
  return template_factory(deserialize_map(), name, args);
}

void EnergyMap::resize_(const int part1_index,
                       const int site1_index,
                       const int part2_index,
                       const int site2_index) {
  FATAL("not implemented");
}
std::vector<double> * EnergyMap::smap_(const int part1_index,
                                       const int site1_index,
                                       const int part2_index,
                                       const int site2_index) {
  FATAL("not implemented");
}
std::vector<double> * EnergyMap::smap_new_(const int part1_index,
                                           const int site1_index,
                                           const int part2_index,
                                           const int site2_index) {
  FATAL("not implemented");
}
const std::vector<std::vector<std::vector<std::vector<std::vector<
    double> > > > >& EnergyMap::map() const {
  FATAL("not implemented");
}

}  // namespace feasst
