#include <cmath>
#include <string>
#include <fstream>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

VisitModelInnerTable::VisitModelInnerTable(argtype * args) : VisitModelInner(args) {
  class_name_ = "VisitModelInnerTable";
  const std::string table_file = str("table_file", args, "-1");
  if (table_file != "-1") {
    read_table_(table_file);
  }
}
VisitModelInnerTable::VisitModelInnerTable(argtype args) : VisitModelInnerTable(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void VisitModelInnerTable::read_table_(const std::string file_name) {
  DEBUG("file_name " << file_name);
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find " << file_name);
  if (file.eof()) {
    return;
  }
  std::string line, descript;
  double double_val;
  int int_val;
  file >> descript >> int_val;
  ASSERT(descript == "site_types", "format error: " << descript);
  const int num_sites = int_val;
  DEBUG("num_sites " << num_sites);

  // size arrays
  site_types_.resize(num_sites);
  resize(num_sites, num_sites, &inner_);
//  resize(num_sites, num_sites, &outer_);
  resize(num_sites, num_sites, &energy_);
  resize(num_sites, num_sites, &delta_);
  resize(num_sites, num_sites, &gamma_);
  for (int type = 0; type < num_sites; ++type) {
    file >> int_val;
    DEBUG("site " << int_val);
    site_types_[type] = int_val;
  }

  for (int itype = 0; itype < num_sites; ++itype) {
    for (int jtype = itype; jtype < num_sites; ++jtype) {
      file >> descript >> int_val;
      ASSERT(descript == "num_orientations_per_pi", "format error: " << descript);
      const int num_orientations_per_pi = int_val;
      DEBUG("num_orientations_per_pi " << num_orientations_per_pi);
      file >> descript >> double_val;
      ASSERT(descript == "gamma", "format error: " << descript);
      const double gamma = double_val;
      DEBUG("gamma " << gamma);
      gamma_[itype][jtype] = gamma;
      file >> descript >> double_val;
      ASSERT(descript == "delta", "format error: " << descript);
      const double delta = double_val;
      DEBUG("delta " << delta);
      delta_[itype][jtype] = delta;
      file >> descript >> int_val;
      ASSERT(descript == "num_z", "format error: " << descript);
      const int num_z = int_val;
      DEBUG("num_z " << num_z);
      int ns1 = 0;
      if (itype == jtype) {
        ns1 = num_orientations_per_pi + 1;
      } else {
        ns1 = 2*num_orientations_per_pi + 1;
      }
      const int ns2 = num_orientations_per_pi + 1;
      const int ne1 = 2*num_orientations_per_pi + 1;
      const int ne2 = ns2;
      const int ne3 = ne1;
      DEBUG("ns1 " << ns1);
      DEBUG("ns2 " << ns2);
      DEBUG("ne1 " << ne1);
      DEBUG("ne2 " << ne2);
      DEBUG("ne3 " << ne3);
      argtype dof = {{"num0", str(ns1)}, {"num1", str(ns2)}, {"num2", str(ne1)},
                     {"num3", str(ne2)}, {"num4", str(ne3)}};
      Table5D * in = &inner_[itype][jtype];
//      Table5D * out = &outer_[itype][jtype];
      Table6D * en = &energy_[itype][jtype];
      *in = Table5D(dof);
      if (num_z > 0) {
//        *out = Table5D(dof);
        dof.insert({"num5", str(num_z)});
        *en = Table6D(dof);
      }
      int num_orientations = 0;
      for (int s1 = 0; s1 < ns1; ++s1) {
      for (int s2 = 0; s2 < ns2; ++s2) {
      for (int e1 = 0; e1 < ne1; ++e1) {
      for (int e2 = 0; e2 < ne2; ++e2) {
      for (int e3 = 0; e3 < ne3; ++e3) {
        ++num_orientations;
        file >> double_val;
        //ASSERT(std::abs(double_val - 1) < NEAR_ZERO, num_orientations << " " << double_val);
        in->set_data(s1, s2, e1, e2, e3, double_val);
        if (num_z > 0) {
          file >> double_val;
          //ASSERT(std::abs(double_val - 1.5) < NEAR_ZERO, num_orientations << " " << double_val);
//          out->set_data(s1, s2, e1, e2, e3, double_val);
          for (int z = 0; z < num_z; ++z) {
            file >> double_val;
            en->set_data(s1, s2, e1, e2, e3, z, double_val);
            //ASSERT(std::abs(double_val) < NEAR_ZERO || std::abs(double_val+1) < NEAR_ZERO,
            //  num_orientations << " " << double_val);
          }
        }
      }}}}}
      INFO("num_orientations " << num_orientations);
    }
  }

  // check for eof
  ASSERT(!file.eof(), "improper table file: " << file_name);
  file >> double_val;
  ASSERT(file.eof(), "improper table file: " << file_name);
}

//bool VisitModelInnerTable::is_outer() const {
//  if (outer_.size() > 0) {
//    if (outer_[0].size() > 0) {
//      if (outer_[0][0].num0() > 1) {
//        return true;
//      }
//    }
//  }
//  return false;
//}

bool VisitModelInnerTable::is_energy_table() const {
  if (energy_.size() > 0) {
    if (energy_[0].size() > 0) {
      if (energy_[0][0].num0() > 1) {
        return true;
      }
    }
  }
  return false;
}

void VisitModelInnerTable::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  aniso_index_ = config->model_params().index("anisotropic");
  DEBUG("aniso_index_ " << aniso_index_);
  for (int t1 = 0; t1 < static_cast<int>(site_types_.size()); ++t1) {
    const int type1 = site_types_[t1];
    for (int t2 = 0; t2 < static_cast<int>(site_types_.size()); ++t2) {
      const int type2 = site_types_[t2];
      const double cutoff = inner_[t1][t2].maximum() + delta_[t1][t2];
      config->set_model_param("cutoff", type1, type2, cutoff);
      INFO("cutoff for " << type1 << "-" << type2 << " site types: " << cutoff);
    }
  }
}

void VisitModelInnerTable::compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) {
  DEBUG("*** VisitModelInnerTable ***");
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  const ModelParam& aniso = model_params.select(aniso_index_);
  int type1 = site1.type();
  int type2 = site2.type();

  // check if sites are anisotropic
  if (aniso.value(type1) < 0.5 || aniso.value(type2) < 0.5) {
    return;
  }
  DEBUG("aniso " << aniso.value(type1) << " " << aniso.value(type2));

  // check if sites are within the global cutoff
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  double squared_distance;
  config->domain().wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance);
  DEBUG("squared_distance " << squared_distance);
  DEBUG("relative " << relative->str());
  if (squared_distance > cutoff*cutoff) {
    return;
  }
  DEBUG("inside global cut");

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
    swap(&type1, &type2);
  }
  DEBUG("flip " << flip);

  // obtain the inverse rotation matrix that sets the reference frame on site 1
  if (flip) {
    site2.euler().compute_rotation_matrix(&rot1_);
  } else {
    site1.euler().compute_rotation_matrix(&rot1_);
  }
  rot1_.transpose();
  DEBUG("rot1 " << rot1_.str());

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
  DEBUG("pos1 " << pos1_.str());
  pos1_.spherical(&sph_);
  DEBUG("sph " << sph_.str());

  // obtain the relative orientation of site 2 in frame of site 1.
  if (flip) {
    site1.euler().compute_rotation_matrix(&rot2_);
  } else {
    site2.euler().compute_rotation_matrix(&rot2_);
  }
  rot1_.multiply(rot2_, &rot3_, &pos1_, &pos2_);
  euler_.set(rot3_);
  DEBUG("euler " << euler_.str());

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
  DEBUG("s1 " << s1 << " s2 " << s2 << " e1 " << e1 << " e2 " << e2 << " e3 " << e3);
  ASSERT(s1 >= 0 && s1 <= 1, "s1: " << s1);
  ASSERT(s2 >= 0 && s2 <= 1, "s2: " << s2);
  ASSERT(e1 >= 0 && e1 <= 1, "e1: " << e1);
  ASSERT(e2 >= 0 && e2 <= 1, "e2: " << e2);
  ASSERT(e3 >= 0 && e3 <= 1, "e3: " << e3);

//  // check the outer cutoff, if applicable.
//  float outer = 0.;
//  const bool global_outer = !is_outer();
//  if (!global_outer) {
//    outer = outer_[type1][type2].linear_interpolation(s1, s2, e1, e2, e3);
//    DEBUG("outer " << outer);
//    if (squared_distance > outer*outer) {
//      return;
//    }
//  }

  // check the inner cutoff.
  const float inner = inner_[type1][type2].linear_interpolation(s1, s2, e1, e2, e3);
  DEBUG("inner " << inner);
  double en;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    DEBUG("hard overlap");
  } else {
    const double gamma = gamma_[type1][type2];
    const double delta = delta_[type1][type2];
    const double outer = inner + delta;
    DEBUG("gamma " << gamma);
    DEBUG("delta " << delta);
    DEBUG("outer " << outer);
    if ((std::abs(gamma) < NEAR_ZERO) && (squared_distance < outer*outer)) {
      en = -1;
    } else if (is_energy_table()) {
      const double rhg = std::pow(inner, gamma);
      const double rcg = std::pow(outer, gamma);
      const double rg = std::pow(squared_distance, 0.5*gamma);
      const double z = std::pow((rg - rcg)/(rhg - rcg), 1./gamma);
      DEBUG("z " << z);
      ASSERT(z >= 0 && z <= 1, "z: " << z);
      en = energy_[type1][type2].linear_interpolation(s1, s2, e1, e2, e3, z);
    } else {
      return;
    }
  }
  DEBUG("en " << en);
  update_ixn(en, part1_index, site1_index, type1, part2_index,
             site2_index, type2, squared_distance, pbc, is_old_config, *config);
}

class MapVisitModelInnerTable {
 public:
  MapVisitModelInnerTable() {
    auto obj = MakeVisitModelInnerTable();
    obj->deserialize_map()["VisitModelInnerTable"] = obj;
  }
};

static MapVisitModelInnerTable mapper_ = MapVisitModelInnerTable();

VisitModelInnerTable::VisitModelInnerTable(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7945, "unrecognized version: " << version);
  feasst_deserialize(&aniso_index_, istr);
}

void VisitModelInnerTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(7945, ostr);
  feasst_serialize(aniso_index_, ostr);
}

double VisitModelInnerTable::second_virial_coefficient(argtype args) const {
  const int expand_t = integer("expand_t", &args, 1);
  const int expand_z = integer("expand_z", &args, 1);
  const int type1 = integer("site_type1", &args, 0);
  const int type2 = integer("site_type2", &args, 0);
  const double beta = dble("beta", &args, 1.);
  FEASST_CHECK_ALL_USED(args);
  if (energy_[0][0].num0() == 1) {
    WARN("Only implemented for hard particles (num_z == 0).");
  }
  double b2 = 0.;
  double b2_h = 0.;
  double b2_a = 0.;
  const Table5D& inner = inner_[type1][type2];
//  const Table5D& outer = outer_[type1][type2];
  const int ns1 = inner.num0();
  const int ns2 = inner.num1();
  const int ne1 = inner.num2();
  const int ne2 = inner.num3();
  const int ne3 = inner.num4();
  DEBUG("ns1 " << ns1 << " ns2 " << ns2 << " ne1 " << ne1 << " ne2 " << ne2 << " ne3 " << ne3);
  const double ds1 = 1./static_cast<double>(ns1*expand_t - 1);
  const double ds2 = 1./static_cast<double>(ns2*expand_t - 1);
  const double de1 = 1./static_cast<double>(ne1*expand_t - 1);
  const double de2 = 1./static_cast<double>(ne2*expand_t - 1);
  const double de3 = 1./static_cast<double>(ne3*expand_t - 1);
  DEBUG("ds1 " << ds1 << " ds2 " << ds2 << " de1 " << de1 << " de2 " << de2 << " de3 " << de3);
  int num_z = 0;
  double dz = 0;
  const Table6D& energy = energy_[type1][type2];
  if (is_energy_table()) {
    num_z = energy.num5();
    dz = 1./static_cast<double>(num_z*expand_z - 1);
  }
  // loop over all orientations
  // loop over corners of multi-dimensional bin for trapezoidal rule (is[1,2], etc)
  #pragma omp parallel for reduction(+:b2_h,b2_a)
  for (int s2 = 0; s2 < ns2*expand_t - 1; ++s2) {
  for (int is2 = 0; is2 <= 1; ++is2) {
  const double sin_s2 = std::sin(PI*(s2+is2)*ds2);
  for (int e2 = 0; e2 < ne2*expand_t - 1; ++e2) {
  for (int ie2 = 0; ie2 <= 1; ++ie2) {
  const double sin_e2 = std::sin(PI*(e2+ie2)*de2);
  for (int s1 = 0; s1 < ns1*expand_t - 1; ++s1) {
  for (int e1 = 0; e1 < ne1*expand_t - 1; ++e1) {
  for (int e3 = 0; e3 < ne3*expand_t - 1; ++e3) {
  for (int is1 = 0; is1 <= 1; ++is1) {
  for (int ie1 = 0; ie1 <= 1; ++ie1) {
  for (int ie3 = 0; ie3 <= 1; ++ie3) {
    // hard particle contribution
    float rh = 0.;
    if (expand_t == 1) {
      rh = inner.data()[s1 + is1][s2 + is2][e1 + ie1][e2 + ie2][e3 + ie3];
    } else {
      rh = inner.linear_interpolation((s1+is1)*ds1, (s2+is2)*ds2, (e1+ie1)*de1, (e2+ie2)*de2, (e3+ie3)*de3);
    }
    b2_h += rh*rh*rh*sin_s2*sin_e2;

    // attractive contribution
    if (num_z > 0) {
      const double rc = rh + delta_[type1][type2];
//      float rc = 0.;
//      if (expand_t == 1) {
//        rc = outer.data()[s1 + is1][s2 + is2][e1 + ie1][e2 + ie2][e3 + ie3];
//      } else {
//        rc = outer.linear_interpolation((s1+is1)*ds1, (s2+is2)*ds2, (e1+ie1)*de1, (e2+ie2)*de2, (e3+ie3)*de3);
//      }
      //ASSERT(std::abs(rc-1.5) < 1e-7, "err");
      for (int z = 0; z < num_z*expand_z - 1; ++z) {
      for (int iz = 0; iz <= 1; ++iz) {
        double u = 0.;
        if (expand_t == 1 && expand_z == 1) {
          u = energy.data()[s1+is1][s2+is2][e1+ie1][e2+ie2][e3+ie3][z+iz];
        } else {
          u = energy.linear_interpolation((s1+is1)*ds1, (s2+is2)*ds2, (e1+ie1)*de1, (e2+ie2)*de2, (e3+ie3)*de3, (z+iz)*dz);
        }
        //ASSERT(std::abs(u+1) < 5e-3, "err");
        const double zval = (z + iz)*dz;
        const double r = (rc - rh)*zval + rh;
        b2_a += (1. - std::exp(-beta*u))*(rc - rh)*r*r*sin_s2*sin_e2;
      }}
    }
  }}}}}}}}}}
  b2_h *= 1./3.; // rh integral prefactor
  b2_a *= dz/2.; // dz norm and extra trapezoid
  b2 += b2_h + b2_a;
  b2 *= ds1*ds2*de1*de2*de3;
  b2 /= std::pow(2., 5); // trapezoid normalization for orientations
  b2 *= 2; // symmetry in s1
  b2 /= 2; // b2 prefactor
  b2 *= PI*PI; // normalization for ds1*ds2
  b2 *= PI/2; // normalization for de1*de2*de3 e.g. 4pi^3/8pi^2
  return b2;
}

}  // namespace feasst
