#include <cmath>  // isnan, pow
#include <string>
#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "server/include/server.h"
#include "server/include/visit_model_inner_server.h"

namespace feasst {

VisitModelInnerServer::VisitModelInnerServer(argtype * args) : VisitModelInner(args) {
  class_name_ = "VisitModelInnerServer";
  server_ = std::make_unique<Server>(args);
  const std::string names = str("server_sites", args, "");
  if (names.empty()) {
    WARN("Deprecate VisitModelInnerServer::server_site[i]->server_sites.");
    int type = 0;
    std::string start = "server_site";
    std::stringstream key;
    key << start << type;
    while (used(key.str(), *args)) {
      site_type_names_.push_back(str(key.str(), args));
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  } else {
    for (const std::string& name : split(names, ',')) {
      site_type_names_.push_back(name);
    }
  }
}
VisitModelInnerServer::VisitModelInnerServer(argtype args) : VisitModelInnerServer(&args) {
  feasst_check_all_used(args);
}
VisitModelInnerServer::~VisitModelInnerServer() {}

void VisitModelInnerServer::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  aniso_index_ = config->model_params().index("anisotropic");
  DEBUG("aniso_index_ " << aniso_index_);
  t2index_.resize(config->num_site_types(), -1);
  for (int t1 = 0; t1 < static_cast<int>(site_type_names_.size()); ++t1) {
    const std::string& stname = site_type_names_[t1];
    const int type1 = config->site_type_name_to_index(stname);
    t2index_[type1] = t1;
  }
}

void VisitModelInnerServer::compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc,
    const double weight) {
  TRACE("*** VisitModelInnerServer ***");
  set_interacted(0);
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  TRACE("aniso_index_ " << aniso_index_);
  const ModelParam& aniso = model_params.select(aniso_index_);
  int type1 = site1.type();
  int type2 = site2.type();

  // check if sites are anisotropic
  TRACE("type1 " << type1)
  TRACE("type2 " << type2);
  TRACE("aniso " << aniso.value(type1) << " " << aniso.value(type2));
  if (aniso.value(type1) < 0.5 || aniso.value(type2) < 0.5) {
    return;
  }

  // check if sites are within the global cutoff
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  double squared_distance;
  config->domain().wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance);
  TRACE("squared_distance " << squared_distance);
  TRACE("relative " << relative->str());
  TRACE("cutoff " << cutoff);
  if (squared_distance > cutoff*cutoff) {
    return;
  }
  TRACE("inside global cut");

  bool flip = false;
  // enforce type1 <= type2 to avoid redundant tables and keep the reference
  // as the smallest possible type.
  // also, if types are the same, and x>0(note,rel is inv), flip 1 and 2
  if (type1 > type2) {
    flip = true;
  } else if (type1 == type2) {
    if (relative->coord(0) > 0) {
      flip = true;
    } else if (relative->coord(0) == 0) {
      if (part1_index > part2_index) {
        flip = true;
      }
    }
  }
  if (flip) {
    feasst_swap(&type1, &type2);
  }
  TRACE("flip " << flip);

  // obtain the inverse rotation matrix that sets the reference frame on site 1
  if (flip) {
    site2.euler().compute_rotation_matrix(&rot1_);
  } else {
    site1.euler().compute_rotation_matrix(&rot1_);
  }
  rot1_.transpose();
  TRACE("rot1 " << rot1_.str());

  // obtain the relative orientation of the centers in spherical coordinates
  if (!flip) {
    relative->multiply(-1); // r12 point toward 1, so reverse direction.
  }
  if (sph_.size() == 0) {
    sph_.set_to_origin(config->dimension());
    pos1_.set_to_origin(config->dimension());
    pos2_.set_to_origin(config->dimension());
  }
  pos1_ = *relative;
  rot1_.multiply(*relative, &pos1_);
  TRACE("pos1 " << pos1_.str());
  pos1_.spherical(&sph_);
  TRACE("sph " << sph_.str());

  // obtain the relative orientation of site 2 in frame of site 1.
  if (flip) {
    site1.euler().compute_rotation_matrix(&rot2_);
  } else {
    site2.euler().compute_rotation_matrix(&rot2_);
  }
  rot1_.multiply(rot2_, &rot3_, &pos1_, &pos2_);
  euler_.set(rot3_);
  TRACE("euler " << euler_.str());

  // obtain the scaled orientational coordinates
  double s1 = 0;
//  s1 = sph_.coord(1)/2/PI;
  if (type1 == type2) {
    s1 = sph_.coord(1)/PI;
  } else {
    s1 = sph_.coord(1)/2/PI;
  }
  if (s1 < 0) {
    s1 += 1;
  }
  double s2 = sph_.coord(2)/PI;
  const double e1 = euler_.phi()/2/PI + 0.5;
  const double e2 = euler_.theta()/PI;
  const double e3 = euler_.psi()/2/PI + 0.5;
  TRACE("s1 " << s1 << " s2 " << s2 << " e1 " << e1 << " e2 " << e2 << " e3 " << e3);
  ASSERT(s1 >= 0 && s1 <= 1, "s1: " << s1);
  ASSERT(s2 >= 0 && s2 <= 1, "s2: " << s2);
  ASSERT(e1 >= 0 && e1 <= 1, "e1: " << e1);
  ASSERT(e2 >= 0 && e2 <= 1, "e2: " << e2);
  ASSERT(e3 >= 0 && e3 <= 1, "e3: " << e3);

  // convert site type to table type
  const int tabtype1 = t2index_[type1];
  const int tabtype2 = t2index_[type2];
  ASSERT(tabtype1 != -1, "site " << type1 << " is anisotropic but not "
    << "included in VisitModelInnerServer.");
  ASSERT(tabtype2 != -1, "site " << type2 << " is anisotropic but not "
    << "included in VisitModelInnerServer.");

  // send the relative orientation to client
  if (!server_->bound()) {
    server_->bind_listen_accept();
  }
  std::stringstream ss;
  ss << squared_distance << "," << s1 << "," << s2 << "," << e1 << "," << e2
     << "," << e3 << "," << type1 << "," << type2;
  // Access neighbor list throught energy_map_
  server_->send(ss.str());
  const int size = server_->receive();
  TRACE(server_->buffer());
  ASSERT(size > 0, "error");
  const double en = weight*std::stod(server_->buffer());
  TRACE("en " << en);
  update_ixn(en, part1_index, site1_index, type1, part2_index,
             site2_index, type2, squared_distance, pbc, is_old_config, *config);
}

FEASST_MAPPER(VisitModelInnerServer, argtype({{"server_sites", "0"}}));

VisitModelInnerServer::VisitModelInnerServer(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 2670 && version <= 2671, "unrecognized version: " << version);
  feasst_deserialize(&aniso_index_, istr);
  if (version >= 2671) {
    feasst_deserialize(&t2index_, istr);
    feasst_deserialize(&site_type_names_, istr);
  }
  feasst_deserialize(server_, istr);
//  feasst_deserialize2(std::move(server_), istr);
  //HWH for unknown reasons, this does not deserialize properly
//  { int existing;
//    istr >> existing;
//    if (existing != 0) {
//      server_ = std::make_unique<Server>(istr);
//    }
//  }
}

void VisitModelInnerServer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(2671, ostr);
  feasst_serialize(aniso_index_, ostr);
  feasst_serialize(t2index_, ostr);
  feasst_serialize(site_type_names_, ostr);
  feasst_serialize(server_, ostr);
}

}  // namespace feasst
