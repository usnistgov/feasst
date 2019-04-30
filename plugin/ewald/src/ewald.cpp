
#include "ewald/include/ewald.h"
#include "core/include/constants.h"

namespace feasst {

void Ewald::set_kmax(const double kmax_squared) {
  kmax_squared_ = kmax_squared;
  kmax_ = static_cast<int>(std::sqrt(kmax_squared_)) + 1;
  kxmax_ = kmax_ + 1;
  kymax_ = 2*kmax_ + 1;
  kzmax_ = kymax_;
//    INFO("kxmax " << kxmax_ << " kymax " << kymax_ << " kzmax " << kzmax_);
}

void Ewald::update_wave_vectors(const Configuration& config) {
  wave_prefactor_.clear();
  wave_num_.clear();
  std::vector<double> kvect(3);
  Position kvec;
  kvec.set_to_origin_3D();
  const Domain& domain = config.domain();
  const double lx = domain.side_length(0);
  const double ly = domain.side_length(1);
  const double lz = domain.side_length(2);
  ASSERT(!domain.is_tilted(), "assumes cuboid domain");
  const double volume = domain.volume();
  const double alpha = config.model_params().property("alpha");
  for (int kx = 0; kx <= kmax_; ++kx) {
  for (int ky = -kmax_; ky <= kmax_; ++ky) {
  for (int kz = -kmax_; kz <= kmax_; ++kz) {
    const int k2 = kx*kx + ky*ky + kz*kz;
    if ( (k2 < kmax_squared_) && (k2 != 0) ) {  // allen tildesley, srsw
    // if ( (k2 <= kmax_squared_) && (k2 != 0) ) {  // gerhard
      kvec.set_vector({2.*PI*kx/lx,
                       2.*PI*ky/ly,
                       2.*PI*kz/lz});
      const double k_sq = kvec.squared_distance();
      double factor = 1.;
      if (kx != 0) factor = 2;
      wave_prefactor_.push_back(2.*PI*factor*exp(-k_sq/4./alpha/alpha)/k_sq/volume);
      //INFO(wave_prefactor_.back() << " k2 " << k_sq << " alpha " << alpha << "  vol " << volume);
      //INFO(2.*PI*factor*exp(-k_sq/4./alpha/alpha));
      wave_num_.push_back(kx);
      wave_num_.push_back(ky);
      wave_num_.push_back(kz);
    }
  }}}
  struct_fact_real_.resize(num_vectors());
  struct_fact_imag_.resize(num_vectors());
}

void Ewald::init_wave_vector_storage(Configuration * config, const int group_index) {
  const Select& selection = config->group_selects()[group_index];
  init_wave_vector_storage(config, selection);
}
void Ewald::init_wave_vector_storage(Configuration * config, const Select& selection) {
  std::stringstream ss;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    for (int site_index : selection.site_indices(select_index)) {
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
            config->add_site_property(ss.str(), 0., part_index, site_index);
          }
        }
      }
    }
  }
}

void Ewald::update_eik(const Select& selection, Configuration * config) {
  const double struct_sign = sign_(selection);
  std::stringstream ss;
  std::string eikrx0("eikrx0");
  const Domain& domain = config->domain();
  const double lx = domain.side_length(0);
  const double ly = domain.side_length(1);
  const double lz = domain.side_length(2);
  const double twopilx = 2.*PI/lx,
               twopily = 2.*PI/ly,
               twopilz = 2.*PI/lz;
  if (selection.trial_state() == "add") {
    init_wave_vector_storage(config, selection);
  }
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    for (int site_index : selection.site_indices(select_index)) {
      // obtain the index for the property
      // this assumes all eik in site are contiguous and ordered
      int eikrx0_index = 0;
      const Site& site = config->select_particle(part_index).site(site_index);
      ASSERT(
        find_in_list(eikrx0, site.properties().property_name(), &eikrx0_index),
        "eikrx0 doesn't exist");
      // calculate eik of kx = 0 explicitly
      ASSERT(kymax_ == kzmax_, "assumption");
      const int eikry0_index = eikrx0_index + kxmax_ + kmax_;
      const int eikrz0_index = eikry0_index + kymax_;
      const int eikix0_index = eikrx0_index + kxmax_ + kymax_ + kzmax_;
      const int eikiy0_index = eikix0_index + kxmax_ + kmax_;
      const int eikiz0_index = eikiy0_index + kymax_;
//        INFO(eikrx0_index << " " << eikry0_index << " " << eikrz0_index << " "
//          << eikix0_index << " " << eikiy0_index << " " << eikiz0_index);
      if (selection.trial_state() != "old") {
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
          const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().property_value();
  //          INFO("test " << eik[eikrx0_index + 1] << " " << cos(twopilx*pos[0]) << " " <<
  //            site.properties().property_value()[0] << " " <<
  //            site.properties().property_value()[eikrx0_index + 1] << " "
  //          );
          config->set_site_property(eikry0_index - 1, eik[eikry0_index + 1], part_index, site_index);
          config->set_site_property(eikiy0_index - 1, -eik[eikiy0_index + 1], part_index, site_index);
          config->set_site_property(eikrz0_index - 1, eik[eikrz0_index + 1], part_index, site_index);
          config->set_site_property(eikiz0_index - 1, -eik[eikiz0_index + 1], part_index, site_index);
        }

        // compute remaining eik by recursion
        for (int kx = 2; kx <= kmax_; ++kx) {
          const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().property_value();
          const double eikr = eik[eikrx0_index + kx - 1]*eik[eikrx0_index + 1] -
            eik[eikix0_index + kx - 1]*eik[eikix0_index + 1];
          config->set_site_property(eikrx0_index + kx, eikr, part_index, site_index);
          const double eiki = eik[eikrx0_index + kx - 1]*eik[eikix0_index + 1] +
            eik[eikix0_index + kx - 1]*eik[eikrx0_index + 1];
          config->set_site_property(eikix0_index + kx, eiki, part_index, site_index);
        }
        for (int ky = 2; ky <= kmax_; ++ky) {
          const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().property_value();
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
          const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().property_value();
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
      const std::vector<double> eik = config->select_particle(part_index).site(site_index).properties().property_value();
      for (int k_index = 0; k_index < num_vectors(); ++k_index) {
        const int kdim = dimension_*k_index;
        const double kx = wave_num_[kdim];
        const double ky = wave_num_[kdim + 1];
        const double kz = wave_num_[kdim + 2];
//          INFO("k " << k_index << " kx " << kx << " ky " << ky << " kz " << kz << " size " << num_vectors() << " kdim " << kdim);
        const double eikrx = eik[eikrx0_index + kx];
        const double eikix = eik[eikix0_index + kx];
        const double eikry = eik[eikry0_index + ky];
        const double eikiy = eik[eikiy0_index + ky];
        const double eikrz = eik[eikrz0_index + kz];
        const double eikiz = eik[eikiz0_index + kz];
//          INFO("eik[r,i]x " << eikrx << " " << eikix << " y " << eikry << " " << eikiy << " z " << eikrz << " " << eikiz << " sz " << eik.size());
        const double eikr = eikrx*eikry*eikrz
                   - eikix*eikiy*eikrz
                   - eikix*eikry*eikiz
                   - eikrx*eikiy*eikiz;
        const double eiki = -eikix*eikiy*eikiz
                   + eikrx*eikry*eikiz
                   + eikrx*eikiy*eikrz
                   + eikix*eikry*eikrz;
//          INFO("charge " << charge << " eikr " << eikr << " eiki " << eiki);
        struct_fact_real_[k_index] += struct_sign*charge*eikr;
        struct_fact_imag_[k_index] += struct_sign*charge*eiki;
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
  feasst_serialize(kmax_, ostr);
  feasst_serialize(kmax_squared_, ostr);
  feasst_serialize(kxmax_, ostr);
  feasst_serialize(kymax_, ostr);
  feasst_serialize(kzmax_, ostr);
  feasst_serialize(wave_prefactor_, ostr);
  feasst_serialize(wave_num_, ostr);
  feasst_serialize(struct_fact_real_, ostr);
  feasst_serialize(struct_fact_imag_, ostr);
  feasst_serialize(struct_fact_real_old_, ostr);
  feasst_serialize(struct_fact_imag_old_, ostr);
  feasst_serialize(stored_energy_, ostr);
  feasst_serialize(stored_energy_old_, ostr);
}

std::shared_ptr<VisitModel> Ewald::create(std::istream& istr) const {
  return std::make_shared<Ewald>(istr);
}

Ewald::Ewald(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(319 == version, version);
  feasst_deserialize(&kmax_, istr);
  feasst_deserialize(&kmax_squared_, istr);
  feasst_deserialize(&kxmax_, istr);
  feasst_deserialize(&kymax_, istr);
  feasst_deserialize(&kzmax_, istr);
  feasst_deserialize(&wave_prefactor_, istr);
  feasst_deserialize(&wave_num_, istr);
  feasst_deserialize(&struct_fact_real_, istr);
  feasst_deserialize(&struct_fact_imag_, istr);
  feasst_deserialize(&struct_fact_real_old_, istr);
  feasst_deserialize(&struct_fact_imag_old_, istr);
  feasst_deserialize(&stored_energy_, istr);
  feasst_deserialize(&stored_energy_old_, istr);
}

}  // namespace feasst
