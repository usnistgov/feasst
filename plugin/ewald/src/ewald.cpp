
#include "ewald/include/ewald.h"
#include "math/include/constants.h"

namespace feasst {

Ewald::Ewald(const argtype& args) {
  class_name_ = "Ewald";
  Arguments args_(args);
  if (args_.key("tolerance").used()) {
    tolerance_ = std::make_shared<double>(args_.dble());
  }
  if (args_.key("alpha").used()) {
    alpha_arg_ = std::make_shared<double>(args_.dble());
  }
  if (args_.key("kxmax").used()) {
    kxmax_arg_ = std::make_shared<int>(args_.integer());
  }
  if (args_.key("kymax").used()) {
    kymax_arg_ = std::make_shared<int>(args_.integer());
  }
  if (args_.key("kzmax").used()) {
    kzmax_arg_ = std::make_shared<int>(args_.integer());
  }
  if (args_.key("kmax_squared").used()) {
    kmax_sq_arg_ = std::make_shared<int>(args_.integer());
  }
}

void Ewald::tolerance_to_alpha_ks(const double tolerance,
    const Configuration& config, double * alpha,
    int * kxmax, int * kymax, int * kzmax, const int num_particles) {
  const double cutoff = config.model_params().cutoff().max();

  // extrapolate to a different number of particles with a factor.
  double factor = 1.;
  if (num_particles != -1) {
    ASSERT(num_particles > 0, "num_particles:" << num_particles
      << " should be > 0");
    factor = num_particles/config.num_particles();
  }

  // determine alpha
  *alpha = std::sqrt(config.num_sites()*cutoff*config.domain()->volume());
  *alpha *= tolerance/(2.*std::sqrt(factor)*sum_squared_charge_(config));
  if (*alpha >= 1.) {
    *alpha = (1.35 - 0.15*log(tolerance))/cutoff; // from LAMMPS
  } else {
    *alpha = std::sqrt(-log(*alpha))/cutoff;
  }
  *kxmax = estimate_kmax_(*alpha, config, tolerance, 0);
  *kymax = estimate_kmax_(*alpha, config, tolerance, 1);
  *kzmax = estimate_kmax_(*alpha, config, tolerance, 2);
}

void Ewald::update_wave_vectors(const Configuration& config) {
  wave_prefactor_.clear();
  wave_num_.clear();
  std::vector<double> kvect(3);
  Position kvec;
  kvec.set_to_origin_3D();
  const Domain * domain = config.domain();
  const double lx = domain->side_length(0);
  const double ly = domain->side_length(1);
  const double lz = domain->side_length(2);
  ASSERT(!domain->is_tilted(), "assumes cuboid domain");
  const double volume = domain->volume();
  const double alpha = config.model_params().property("alpha");
  for (int kx = 0; kx <= kmax_; ++kx) {
  for (int ky = -kmax_; ky <= kmax_; ++ky) {
  for (int kz = -kmax_; kz <= kmax_; ++kz) {
    kvec.set_vector({2.*PI*kx/lx,
                     2.*PI*ky/ly,
                     2.*PI*kz/lz});
    const double k_sq = kvec.squared_distance();
    if ( (k_sq < kmax_squared_) and (std::abs(k_sq) > NEAR_ZERO) ) { // allen tildesley, srsw
    // if ( (k2 <= kmax_squared_) and (std::abs(k_sq) > NEAR_ZERO) ) {  // gerhard
      double factor = 1.;
      if (kx != 0) factor = 2;
      wave_prefactor_.push_back(2.*PI*factor*exp(-k_sq/4./alpha/alpha)/k_sq/volume);
      DEBUG(wave_prefactor_.back() << " k2 " << k_sq << " alpha " << alpha << "  vol " << volume);
      DEBUG(2.*PI*factor*exp(-k_sq/4./alpha/alpha));
      wave_num_.push_back(kx);
      wave_num_.push_back(ky);
      wave_num_.push_back(kz);
    }
  }}}
  struct_fact_real_.resize(num_vectors());
  struct_fact_imag_.resize(num_vectors());
  struct_fact_real_new_.resize(num_vectors());
  struct_fact_imag_new_.resize(num_vectors());
}

std::vector<std::string> Ewald::eik_gen_() {
  std::vector<std::string> eiks;
  std::stringstream ss;
  for (std::string comp : {"r", "i"}) {
    for (std::string coord : {"x", "y", "z"}) {
      int num_k = -1;
      if (coord == "x") {
        num_k = kxmax_;
      } else if (coord == "y") {
        num_k = kymax_;
      } else if (coord == "z") {
        num_k = kzmax_;
      } else {
        ERROR("unrecognized coord: " << coord);
      }
      for (int k = 0; k < num_k; ++k) {
        ss.str("");
        ss << "eik" << comp << coord << k;
        eiks.push_back(ss.str());
      }
    }
  }
  return eiks;
}

void Ewald::init_wave_vector_storage(Configuration * config, const int group_index) {
  const Select& selection = config->group_selects()[group_index];
  init_wave_vector_storage(config, selection);
  // also add eik properties to the particle types
  const std::vector<std::string> eiks = eik_gen_();
  for (int part_type = 0; part_type < config->num_particle_types();
       ++part_type) {
    const Particle& part = config->particle_type(part_type);
    for (int site = 0; site < part.num_sites(); ++site) {
      for (const std::string& eik : eiks) {
        config->add_or_set_particle_type_site_property(eik, 0.,
          part_type, site);
      }
    }
  }
}

void Ewald::init_wave_vector_storage(Configuration * config, const Select& selection) {
  const std::vector<std::string> eiks = eik_gen_();
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    for (int site_index : selection.site_indices(select_index)) {
      for (const std::string& eik : eiks) {
        config->add_or_set_site_property(eik, 0., part_index, site_index);
      }
    }
  }
}

void Ewald::precompute(Configuration * config) {
  if (tolerance_) {
    ASSERT(!alpha_arg_ && !kxmax_arg_ && !kymax_arg_ && !kzmax_arg_,
      "tolerance overrides all other arguments");
    double alpha;
    tolerance_to_alpha_ks(*tolerance_, *config, &alpha, &kxmax_, &kymax_, &kzmax_);
    config->add_or_set_model_param("alpha", alpha);
  } else if (kmax_sq_arg_) {
    ASSERT(!kxmax_arg_ && !kymax_arg_ && !kzmax_arg_,
      "kmax_squared argument overrides k[x,y,z]max arguments.");
    ASSERT(config->domain()->is_cubic(),
      "Domain must be cubic to set kmax_squared");
    kmax_ = static_cast<int>(std::sqrt(*kmax_sq_arg_)) + 1;
    kxmax_ = kmax_ + 1;
    kymax_ = 2*kmax_ + 1;
    kzmax_ = kymax_;
    kmax_squared_ = *kmax_sq_arg_*pow(2.*PI/config->domain()->min_side_length(), 2);
    DEBUG("kxmax " << kxmax_ << " kymax " << kymax_ << " kzmax " << kzmax_);
    DEBUG("kmax_squared_ " << kmax_squared_);
    ASSERT(alpha_arg_, "alpha is required with kmax_squared argument");
    config->add_or_set_model_param("alpha", *alpha_arg_);
  } else {
    kxmax_ = *kxmax_arg_;
    kymax_ = *kymax_arg_;
    kzmax_ = *kzmax_arg_;
    double gsqxmx = pow(2*PI*kxmax_/config->domain()->side_length(0), 2);
    double gsqymx = pow(2*PI*kymax_/config->domain()->side_length(1), 2);
    double gsqzmx = pow(2*PI*kzmax_/config->domain()->side_length(2), 2);
    kmax_squared_ = std::max(gsqxmx, gsqymx);
    kmax_squared_ = std::max(kmax_squared_, gsqzmx);
    ASSERT(alpha_arg_ && kxmax_arg_ && kymax_arg_ && kzmax_arg_,
      "if tolerance or kmax_squared are not given, then alpha, kxmax, kymax "
      << "and kzmax are required");
    config->add_or_set_model_param("alpha", *alpha_arg_);
  }
  config->add_excluded_property("eik");
  update_wave_vectors(*config);
  init_wave_vector_storage(config);
}

void Ewald::update_struct_fact_eik(const Select& selection,
                                   Configuration * config,
                                   std::vector<double> * struct_fact_real,
                                   std::vector<double> * struct_fact_imag) {
  ASSERT(struct_fact_real->size() == struct_fact_real_.size(),
    "While struct_fact_real_ is of size: " << struct_fact_real_.size() <<
    " struct_fact_real is of size: " << struct_fact_real->size());
  std::stringstream ss;
  const Domain * domain = config->domain();
  const double lx = domain->side_length(0);
  const double ly = domain->side_length(1);
  const double lz = domain->side_length(2);
  const double twopilx = 2.*PI/lx,
               twopily = 2.*PI/ly,
               twopilz = 2.*PI/lz;
  const int state = selection.trial_state();
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const double struct_sign = sign_(selection, select_index);
    for (int ss_index = 0; ss_index < selection.num_sites(select_index); ++ss_index) {
      const int site_index = selection.site_index(select_index, ss_index);
      // obtain the index for the property
      // this assumes all eik in site are contiguous and ordered
      const Site& site = config->select_particle(part_index).site(site_index);
      if (site.is_physical()) {
        const int eikrx0_index = find_eikrx0_(site);
        // calculate eik of kx = 0 explicitly
        ASSERT(kymax_ == kzmax_, "assumption");
        const int eikry0_index = eikrx0_index + kxmax_ + kmax_;
        const int eikrz0_index = eikry0_index + kymax_;
        const int eikix0_index = eikrx0_index + kxmax_ + kymax_ + kzmax_;
        const int eikiy0_index = eikix0_index + kxmax_ + kmax_;
        const int eikiz0_index = eikiy0_index + kymax_;
        TRACE(eikrx0_index << " " << eikry0_index << " " << eikrz0_index << " "
          << eikix0_index << " " << eikiy0_index << " " << eikiz0_index);

        // update the eik of the selection
        if (state != 0 && state != 2) {
          config->set_site_property(eikrx0_index, 1., part_index, site_index);
          config->set_site_property(eikix0_index, 0., part_index, site_index);
          config->set_site_property(eikry0_index, 1., part_index, site_index);
          config->set_site_property(eikiy0_index, 0., part_index, site_index);
          config->set_site_property(eikrz0_index, 1., part_index, site_index);
          config->set_site_property(eikiz0_index, 0., part_index, site_index);

          // calculate eik of kx = +/-1 explicitly
          const std::vector<double> pos = config->select_particle(part_index).site(site_index).position().coord();
          config->set_site_property(eikrx0_index + 1, cos(twopilx*pos[0]), part_index, site_index);
          config->set_site_property(eikix0_index + 1, sin(twopilx*pos[0]), part_index, site_index);
          config->set_site_property(eikry0_index + 1, cos(twopily*pos[1]), part_index, site_index);
          config->set_site_property(eikiy0_index + 1, sin(twopily*pos[1]), part_index, site_index);
          config->set_site_property(eikrz0_index + 1, cos(twopilz*pos[2]), part_index, site_index);
          config->set_site_property(eikiz0_index + 1, sin(twopilz*pos[2]), part_index, site_index);
          {
            const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().values();
            TRACE("test " << eik[eikrx0_index + 1] << " " << cos(twopilx*pos[0]) << " " <<
              site.properties().values()[0] << " " <<
              site.properties().values()[eikrx0_index + 1] << " "
            );
            config->set_site_property(eikry0_index - 1, eik[eikry0_index + 1], part_index, site_index);
            config->set_site_property(eikiy0_index - 1, -eik[eikiy0_index + 1], part_index, site_index);
            config->set_site_property(eikrz0_index - 1, eik[eikrz0_index + 1], part_index, site_index);
            config->set_site_property(eikiz0_index - 1, -eik[eikiz0_index + 1], part_index, site_index);
          }

          // compute remaining eik by recursion
          for (int kx = 2; kx <= kmax_; ++kx) {
            const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().values();
            const double eikr = eik[eikrx0_index + kx - 1]*eik[eikrx0_index + 1] -
              eik[eikix0_index + kx - 1]*eik[eikix0_index + 1];
            config->set_site_property(eikrx0_index + kx, eikr, part_index, site_index);
            const double eiki = eik[eikrx0_index + kx - 1]*eik[eikix0_index + 1] +
              eik[eikix0_index + kx - 1]*eik[eikrx0_index + 1];
            config->set_site_property(eikix0_index + kx, eiki, part_index, site_index);
          }
          for (int ky = 2; ky <= kmax_; ++ky) {
            const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().values();
            const double eikr = eik[eikry0_index + ky - 1]*eik[eikry0_index + 1] -
              eik[eikiy0_index + ky - 1]*eik[eikiy0_index + 1];
            config->set_site_property(eikry0_index + ky, eikr, part_index, site_index);
            const double eiki = eik[eikry0_index + ky - 1]*eik[eikiy0_index + 1] +
              eik[eikiy0_index + ky - 1]*eik[eikry0_index + 1];
            config->set_site_property(eikiy0_index + ky, eiki, part_index, site_index);
            config->set_site_property(eikry0_index - ky, eikr, part_index, site_index);
            config->set_site_property(eikiy0_index - ky, -eiki, part_index, site_index);
          }
          for (int kz = 2; kz <= kmax_; ++kz) {
            const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().values();
            const double eikr = eik[eikrz0_index + kz - 1]*eik[eikrz0_index + 1] -
              eik[eikiz0_index + kz - 1]*eik[eikiz0_index + 1];
            config->set_site_property(eikrz0_index + kz, eikr, part_index, site_index);
            const double eiki = eik[eikrz0_index + kz - 1]*eik[eikiz0_index + 1] +
              eik[eikiz0_index + kz - 1]*eik[eikrz0_index + 1];
            config->set_site_property(eikiz0_index + kz, eiki, part_index, site_index);
            config->set_site_property(eikrz0_index - kz, eikr, part_index, site_index);
            config->set_site_property(eikiz0_index - kz, -eiki, part_index, site_index);
          }
        }

        // compute structure factor
        const int type = site.type();
        const double charge = config->model_params().charge().value(type);
        const Particle& p = config->select_particle(part_index);
        const std::vector<double>& eik = p.site(site_index).properties().values();
        for (int k_index = 0; k_index < num_vectors(); ++k_index) {
          const int kdim = dimension_*k_index;
          const double kx = wave_num_[kdim];
          const double ky = wave_num_[kdim + 1];
          const double kz = wave_num_[kdim + 2];
          TRACE("k " << k_index << " kx " << kx << " ky " << ky << " kz " << kz << " size " << num_vectors() << " kdim " << kdim);
          const double eikrx = eik[eikrx0_index + kx];
          const double eikix = eik[eikix0_index + kx];
          const double eikry = eik[eikry0_index + ky];
          const double eikiy = eik[eikiy0_index + ky];
          const double eikrz = eik[eikrz0_index + kz];
          const double eikiz = eik[eikiz0_index + kz];
          TRACE("eik[r,i]x " << eikrx << " " << eikix << " y " << eikry << " " << eikiy << " z " << eikrz << " " << eikiz << " sz " << eik.size());
          const double eikr = eikrx*eikry*eikrz
                     - eikix*eikiy*eikrz
                     - eikix*eikry*eikiz
                     - eikrx*eikiy*eikiz;
          const double eiki = -eikix*eikiy*eikiz
                     + eikrx*eikry*eikiz
                     + eikrx*eikiy*eikrz
                     + eikix*eikry*eikrz;
          TRACE("charge " << charge << " eikr " << eikr << " eiki " << eiki << " sign " << struct_sign);
          TRACE(struct_fact_real->size());
          (*struct_fact_real)[k_index] += struct_sign*charge*eikr;
          (*struct_fact_imag)[k_index] += struct_sign*charge*eiki;
        }
      }
    }
  }
}

class MapEwald {
 public:
  MapEwald() {
    Ewald().deserialize_map()["Ewald"] = std::make_shared<Ewald>();
  }
};

static MapEwald mapper_ = MapEwald();

void Ewald::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(319, ostr);
  feasst_serialize_sp(tolerance_, ostr);
  feasst_serialize_sp(alpha_arg_, ostr);
  feasst_serialize_sp(kxmax_arg_, ostr);
  feasst_serialize_sp(kymax_arg_, ostr);
  feasst_serialize_sp(kzmax_arg_, ostr);
  feasst_serialize_sp(kmax_sq_arg_, ostr);
  feasst_serialize(kmax_, ostr);
  feasst_serialize(kmax_squared_, ostr);
  feasst_serialize(kxmax_, ostr);
  feasst_serialize(kymax_, ostr);
  feasst_serialize(kzmax_, ostr);
  feasst_serialize(wave_prefactor_, ostr);
  feasst_serialize(wave_num_, ostr);
  feasst_serialize(struct_fact_real_, ostr);
  feasst_serialize(struct_fact_imag_, ostr);
  feasst_serialize(struct_fact_real_new_, ostr);
  feasst_serialize(struct_fact_imag_new_, ostr);
  feasst_serialize(stored_energy_, ostr);
}

std::shared_ptr<VisitModel> Ewald::create(std::istream& istr) const {
  return std::make_shared<Ewald>(istr);
}

Ewald::Ewald(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(319 == version, version);
//  feasst_deserialize(tolerance_, istr);
//  feasst_deserialize(alpha_arg_, istr);
  double value;
  int existing;
  istr >> existing;
  if (existing != 0) {
    istr >> value;
    tolerance_ = std::make_shared<double>(value);
  }
  istr >> existing;
  if (existing != 0) {
    istr >> value;
    alpha_arg_ = std::make_shared<double>(value);
  }
  //feasst_deserialize(kxmax_arg_, istr);
  //feasst_deserialize(kymax_arg_, istr);
  //feasst_deserialize(kzmax_arg_, istr);
  //feasst_deserialize(kmax_sq_arg_, istr);
  int int_value;
  istr >> existing;
  if (existing != 0) {
    istr >> int_value;
    kxmax_arg_ = std::make_shared<int>(int_value);
  }
  istr >> existing;
  if (existing != 0) {
    istr >> int_value;
    kymax_arg_ = std::make_shared<int>(int_value);
  }
  istr >> existing;
  if (existing != 0) {
    istr >> int_value;
    kzmax_arg_ = std::make_shared<int>(int_value);
  }
  istr >> existing;
  if (existing != 0) {
    istr >> int_value;
    kmax_sq_arg_ = std::make_shared<int>(int_value);
  }
  feasst_deserialize(&kmax_, istr);
  feasst_deserialize(&kmax_squared_, istr);
  feasst_deserialize(&kxmax_, istr);
  feasst_deserialize(&kymax_, istr);
  feasst_deserialize(&kzmax_, istr);
  feasst_deserialize(&wave_prefactor_, istr);
  feasst_deserialize(&wave_num_, istr);
  feasst_deserialize(&struct_fact_real_, istr);
  feasst_deserialize(&struct_fact_imag_, istr);
  feasst_deserialize(&struct_fact_real_new_, istr);
  feasst_deserialize(&struct_fact_imag_new_, istr);
  feasst_deserialize(&stored_energy_, istr);
}

}  // namespace feasst
