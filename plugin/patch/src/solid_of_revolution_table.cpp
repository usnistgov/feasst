#include <cmath>  // isnan, pow
#include <string>
#include <fstream>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "patch/include/solid_of_revolution_table.h"

namespace feasst {

SolidOfRevolutionTable::SolidOfRevolutionTable(argtype * args) : VisitModelInner(args) {
  class_name_ = "SolidOfRevolutionTable";
  table_file_ = str("table_file", args, "");
  ignore_energy_ = boolean("ignore_energy", args, false);
}
SolidOfRevolutionTable::SolidOfRevolutionTable(argtype args) : SolidOfRevolutionTable(&args) {
  feasst_check_all_used(args);
}

void SolidOfRevolutionTable::read_table_(const std::string file_name,
    const bool ignore_energy,
    Configuration * config) {
  TRACE("file_name " << file_name);
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find table file:" << file_name);
  if (file.eof()) {
    return;
  }
  std::string line, descript;
  double double_val;
  int int_val;
  file >> descript >> int_val;
  ASSERT(descript == "site_types", "format error: " << descript);
  const int num_sites = int_val;
  TRACE("num_sites " << num_sites);

  // size arrays
  site_types_.resize(num_sites);
  std::vector<std::vector<std::shared_ptr<Table3D> > > * inner = config->get_table3d();
  std::vector<std::vector<std::shared_ptr<Table4D> > > * energy = config->get_table4d();
  resize(num_sites, num_sites, inner);
  resize(num_sites, num_sites, energy);
  for (int i = 0; i < num_sites; ++i) {
    for (int j = 0; j < num_sites; ++j) {
      (*inner)[i][j] = std::make_shared<Table3D>();
      (*energy)[i][j] = std::make_shared<Table4D>();
    }
  }
  resize(num_sites, num_sites, &delta_);
  resize(num_sites, num_sites, &gamma_);
  resize(num_sites, num_sites, &smoothing_distance_);
  for (int type = 0; type < num_sites; ++type) {
    file >> int_val;
    TRACE("site " << int_val);
    site_types_[type] = int_val;
  }

  for (int itype = 0; itype < num_sites; ++itype) {
    for (int jtype = itype; jtype < num_sites; ++jtype) {
      file >> descript >> int_val;
      ASSERT(descript == "num_orientations_per_half_pi", "format error: " << descript);
      const int num_orientations_per_half_pi = int_val;
      TRACE("num_orientations_per_half_pi " << num_orientations_per_half_pi);
      file >> descript >> double_val;
      ASSERT(descript == "gamma", "format error: " << descript);
      const double gamma = double_val;
      TRACE("gamma " << gamma);
      gamma_[itype][jtype] = gamma;
      file >> descript >> double_val;
      ASSERT(descript == "delta", "format error: " << descript);
      const double delta = double_val;
      TRACE("delta " << delta);
      delta_[itype][jtype] = delta;
      file >> descript >> int_val;
      ASSERT(descript == "num_z", "format error: " << descript);
      const int num_z = int_val;
      TRACE("num_z " << num_z);
      file >> descript >> double_val;
      ASSERT(descript == "smoothing_distance", "format error: " << descript);
      smoothing_distance_[itype][jtype] = double_val;
      TRACE("smoothing_distance " << smoothing_distance_[itype][jtype]);
      const int nt1 = num_orientations_per_half_pi + 1;
      const int nt2 = num_orientations_per_half_pi + 1;
      const int np = 2*num_orientations_per_half_pi + 1;
      TRACE("nt1 " << nt1);
      TRACE("nt2 " << nt2);
      TRACE("np " << np);
      argtype dof = {{"num0", str(nt1)}, {"num1", str(nt2)}, {"num2", str(np)}};
      Table3D * in = (*inner)[itype][jtype].get();
      Table4D * en = (*energy)[itype][jtype].get();
      *in = Table3D(dof);
      if (num_z > 0 && !ignore_energy) {
        dof.insert({"num3", str(num_z)});
        *en = Table4D(dof);
      }
      int num_orientations = 0;
      for (int t1 = 0; t1 < nt1; ++t1) {
      for (int t2 = 0; t2 < nt2; ++t2) {
      for (int p = 0; p < np; ++p) {
        ++num_orientations;
      }}}
      TRACE("num_orientations " << num_orientations);
      // temporarily store values to recall for redundant orientations
      std::vector<std::vector<double> > tmp_store;
      resize(num_orientations, num_z + 2, &tmp_store);
      int ior = 0;
      for (int t1 = 0; t1 < nt1; ++t1) {
      for (int t2 = 0; t2 < nt2; ++t2) {
      for (int p = 0; p < np; ++p) {
        file >> double_val;
        // check for unique orientations
        if (std::abs(double_val + 1) < NEAR_ZERO) {
          int unique_ior;
          file >> unique_ior;
          in->set_data(t1, t2, p, tmp_store[unique_ior][0]);
          if (num_z > 0) {
            for (int z = 0; z < num_z; ++z) {
              if (!ignore_energy) {
                en->set_data(t1, t2, p, z, tmp_store[unique_ior][z + 1]);
              }
            }
          }
        } else {
          in->set_data(t1, t2, p, double_val);
          tmp_store[ior][0] = double_val;
          if (num_z > 0) {
            for (int z = 0; z < num_z; ++z) {
              file >> double_val;
              if (!ignore_energy) {
                en->set_data(t1, t2, p, z, double_val);
                tmp_store[ior][z + 1] = double_val;
              }
            }
          }
        }
        ++ior;
      }}}

      // check the table for bad values
      ASSERT(!has_bad_value(in->data()), "error");
      ASSERT(!has_bad_value(en->data()), "error");
    }
  }

  // check for eof
  ASSERT(!file.eof(), "improper table file: " << file_name);
  file >> double_val;
  ASSERT(file.eof(), "improper table file: " << file_name);
}

bool SolidOfRevolutionTable::is_energy_table(const std::vector<std::vector<std::shared_ptr<Table4D> > >& energy) const {
  if (energy.size() > 0) {
    if (energy[0].size() > 0) {
      if (energy[0][0]->num0() > 1) {
        return true;
      }
    }
  }
  return false;
}

void SolidOfRevolutionTable::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  read_table_(table_file_, ignore_energy_, config);
  director_index_ = config->model_params().index("director");
  TRACE("director_index_ " << director_index_);
  t2index_.resize(config->num_site_types(), 0);
  const std::vector<std::vector<std::shared_ptr<Table3D> > >& inner = config->table3d();
  for (int t1 = 0; t1 < static_cast<int>(site_types_.size()); ++t1) {
    const int type1 = site_types_[t1];
    ASSERT(type1 < config->num_site_types(),"site type: " << type1 <<
      " in table > number of site types:" << config->num_site_types());
    t2index_[type1] = t1;
    for (int t2 = t1; t2 < static_cast<int>(site_types_.size()); ++t2) {
      const int type2 = site_types_[t2];
      ASSERT(type2 < config->num_site_types(),"site type: " << type2 <<
        " in table > number of site types:" << config->num_site_types());
      const double cutoff = inner[t1][t2]->maximum() + delta_[t1][t2];
      // HWH this assumes director types are always following a center type
      config->set_model_param("cutoff", type1-1, type2-1, cutoff);
      config->set_model_param("cutoff", type2-1, type1-1, cutoff);
      config->set_model_param("cutoff", type1, type2, cutoff);
      config->set_model_param("cutoff", type2, type1, cutoff);
      INFO("cutoff for " << type1-1 << "-" << type2-1 << " site types: " << cutoff);
      INFO("cutoff for " << type2-1 << "-" << type1-1 << " site types: " << cutoff);
      INFO("cutoff for " << type1 << "-" << type2 << " site types: " << cutoff);
      INFO("cutoff for " << type2 << "-" << type1 << " site types: " << cutoff);
    }
  }
  dir1_pos_.set_to_origin(config->domain().dimension());
  dir2_pos_.set_to_origin(config->domain().dimension());
}

void SolidOfRevolutionTable::compute(
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
  TRACE("*** SolidOfRevolutionTable ***");
  set_interacted(0);
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  TRACE("site1.pos " << site1.position().str());
  TRACE("site2.pos " << site2.position().str());
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  TRACE("director_index_ " << director_index_);
  const ModelParam& director = model_params.select(director_index_);
  int type1 = site1.type();
  int type2 = site2.type();

  double squared_distance;
  config->domain().wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance);
  TRACE("squared_distance " << squared_distance);
  TRACE("relative " << relative->str());

  // loop through bonds on center1
  const int part1_type = part1.type();
  for (const int dir1_index : config->particle_type(part1_type).bond_neighbors(site1_index)) {
    TRACE("dir1_index " << dir1_index);
    const Site& dir1 = part1.site(dir1_index);
    if (dir1.is_physical()) {
      const int dir1_type = dir1.type();
      if (director.value(dir1_type) > 0.5) {
        TRACE("dir1 pos " << dir1.position().str());
        // loop through bonds on center2
        const int part2_type = part2.type();
        for (const int dir2_index : config->particle_type(part2_type).bond_neighbors(site2_index)) {
          TRACE("dir2_index " << dir2_index);
          const Site& dir2 = part2.site(dir2_index);
          if (dir2.is_physical()) {
            const int dir2_type = dir2.type();
            TRACE("dir2_type " << dir2_type);
            if (director.value(dir2_type) > 0.5) {
              TRACE("dir2 pos " << dir2.position().str());
              const double cutoff = model_params.select(cutoff_index()).mixed_values()[dir1_type][dir2_type];
              TRACE("cutoff " << cutoff);
              if (squared_distance < cutoff*cutoff) {
                dir1_pos_.set_vector(dir1.position().coord());
                dir1_pos_.subtract(site1.position());
                dir1_pos_.normalize();
                TRACE("dir1_pos " << dir1_pos_.str());
                dir2_pos_.set_vector(dir2.position().coord());
                dir2_pos_.subtract(site2.position());
                dir2_pos_.normalize();
                TRACE("dir2_pos " << dir2_pos_.str());

                // for i!=j, there is no i-j swap symmetry
                // tables are referenced such that i is never greater than j
                // if this is the case, then swap i and j
                int type1tmp = type1;
                int type2tmp = type2;
                int site1_indextmp = site1_index;
                int site2_indextmp = site2_index;
                int part1_indextmp = part1_index;
                int part2_indextmp = part2_index;
                if (type1 > type2) {
                  relative->multiply(-1);
                  type1tmp = type2;
                  type2tmp = type1;
                  site1_indextmp = site2_index;
                  site2_indextmp = site1_index;
                  part1_indextmp = part2_index;
                  part2_indextmp = part1_index;
                }

                // uij is unit vector of xij
                const double r = std::sqrt(squared_distance);
                const double ux = relative->coord(0)/r;
                const double uy = relative->coord(1)/r;
                const double uz = relative->coord(2)/r;
                TRACE("ux " << ux << " uy " << uy << " uz " << uz);

                // cosine of thetai, angle between uij and imol
                double i2x = dir1_pos_.coord(0);
                double i2y = dir1_pos_.coord(1);
                double i2z = dir1_pos_.coord(2);
                double cosi = -(i2x*ux+i2y*uy+i2z*uz);
                TRACE("cosi " << cosi);

                // compute symmetry operations on thetai
                //  alternatively, if (cosi<0), cosi = -cosi
                if (cosi < 0) {
                  cosi  = -cosi;
                  i2x = -i2x;
                  i2y = -i2y;
                  i2z = -i2z;
                }
                // note that its possible cosi could be == 1+1e-15. woudl this segfault?

                // cosine of thetaj, angle between uij and jmol
                double j2x = dir2_pos_.coord(0);
                double j2y = dir2_pos_.coord(1);
                double j2z = dir2_pos_.coord(2);
                TRACE("j2x " << j2x);
                TRACE("j2y " << j2y);
                TRACE("j2z " << j2z);
                double cosj = j2x*ux+j2y*uy+j2z*uz;
                TRACE("cosj " << cosj);

                // compute symmetry operations on thetai
                if (cosj < 0) {
                  cosj  = -cosj;
                  j2x = -j2x;
                  j2y = -j2y;
                  j2z = -j2z;
                }

                // ni, normal of imol and uij
                const double nix = i2y*uz - i2z*uy,
                             niy = i2z*ux - i2x*uz,
                             niz = i2x*uy - i2y*ux,
                  nis = std::sqrt(nix*nix+niy*niy+niz*niz);
                double nixn = nix/nis,
                  niyn = niy/nis,
                  nizn = niz/nis;
                if (nis == 0) {
                  nixn = 1.;
                  niyn = nizn = 0.;
                }

                // nj, normal of jmol and -uij
                const double njx = -(j2y*uz - j2z*uy),
                             njy = -(j2z*ux - j2x*uz),
                             njz = -(j2x*uy - j2y*ux),
                  njs = std::sqrt(njx*njx+njy*njy+njz*njz);
                double njxn = njx/njs,
                  njyn = njy/njs,
                  njzn = njz/njs;
                if (njs == 0) {
                  njxn = 1.;
                  njyn = njzn = 0.;
                }

                // cosine of psi is the dot product of ni and nj
                //  no symmetry operation on psi necessary
                const double cospsi = nixn*njxn+niyn*njyn+nizn*njzn;
                TRACE("cospsi " << cospsi);
                const double half_cospsi_p_one = 0.5*(cospsi+1);
                TRACE("half_cospsi_p_one " << half_cospsi_p_one);

                // use i-j swap symmetry to always access the same table element
                //  whether or its i-j or j-i
                //  because otherwise, the monte carlo volume integration could
                //  lead to noise in the potential
                if ( (type1 == type2) && (cosi > cosj) ) {
                  const double cositmp = cosi;
                  cosi = cosj;
                  cosj = cositmp;
                }

                TRACE("cosi " << cosi << " cosj " << cosj << " half_cospsi_p_one " << half_cospsi_p_one);

                // convert site type to table type
                const int tabtype1 = t2index_[type1tmp];
                const int tabtype2 = t2index_[type2tmp];

                // check the inner cutoff.
                const std::vector<std::vector<std::shared_ptr<Table3D> > >& innert = config->table3d();
                TRACE("size1 " << innert.size());
                TRACE("size2 " << innert[0].size());
                double en = 0.;
                if (squared_distance < 1e-6) {
                  en = NEAR_INFINITY;
                  TRACE("hard overlap");
                } else {
                  const float inner = innert[tabtype1][tabtype2]->linear_interpolation(cosi, cosj, half_cospsi_p_one);
                  TRACE("inner " << inner);
                  if (squared_distance < inner*inner) {
                    en = NEAR_INFINITY;
                    TRACE("hard overlap");
                  } else {
                    const double delta = delta_[tabtype1][tabtype2];
                    const double outer = inner + delta;
                    TRACE("delta " << delta);
                    TRACE("outer " << outer);
                    if (squared_distance < outer*outer) {
                      const double gamma = gamma_[tabtype1][tabtype2];
                      TRACE("gamma " << gamma);
                      const std::vector<std::vector<std::shared_ptr<Table4D> > >& energyt = config->table4d();
                      if ((std::abs(gamma) < NEAR_ZERO)) {
                        en = -1;
                      } else if (is_energy_table(energyt)) {
                        const double smooth = smoothing_distance_[tabtype1][tabtype2];
                        const double rhg = std::pow(inner, gamma);
                        const double rcg = std::pow(outer - smooth, gamma);
                        const double rg = std::pow(squared_distance, 0.5*gamma);
                        double z = (rg - rhg)/(rcg - rhg);
                        if (z < 0 && z > -1e-6) {
                          z = 0.;
                        }
                        TRACE("z " << z);
                        if (z > 1.) {
                          en = energyt[tabtype1][tabtype2]->linear_interpolation(cosi, cosj, half_cospsi_p_one, 1.);
                          const double dx = outer - std::sqrt(squared_distance);
                          TRACE("dx " << dx);
                          if (dx > smooth && dx < smooth + 1e-5) {
                            en = 0.;
                          } else {
                            ASSERT(dx >= 0 && dx <= smooth, "dx: " << MAX_PRECISION << dx);
                            en *= dx/smooth;
                          }
                        } else {
                          ASSERT(z >= 0 && z <= 1, "z: " << MAX_PRECISION << z);
                          en = energyt[tabtype1][tabtype2]->linear_interpolation(cosi, cosj, half_cospsi_p_one, z);
                        }
                      } else {
                        return;
                      }
                    }
                  }
                  en *= weight;
                  TRACE("en " << en);
                  update_ixn(en, part1_indextmp, site1_indextmp, type1tmp, part2_indextmp,
                             site2_indextmp, type2tmp, squared_distance, pbc, is_old_config, *config);
                }
              }
            }
          }
        }
      }
    }
  }
}

FEASST_MAPPER(SolidOfRevolutionTable,);

SolidOfRevolutionTable::SolidOfRevolutionTable(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2304, "unrecognized version: " << version);
  feasst_deserialize(&director_index_, istr);
  feasst_deserialize(&table_file_, istr);
  feasst_deserialize(&ignore_energy_, istr);
  feasst_deserialize_fstobj(&dir1_pos_, istr);
  feasst_deserialize_fstobj(&dir2_pos_, istr);
}

void SolidOfRevolutionTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(2304, ostr);
  feasst_serialize(director_index_, ostr);
  feasst_serialize(table_file_, ostr);
  feasst_serialize(ignore_energy_, ostr);
  feasst_serialize_fstobj(dir1_pos_, ostr);
  feasst_serialize_fstobj(dir2_pos_, ostr);
}

}  // namespace feasst
