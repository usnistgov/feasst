#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/accumulator.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/scattering.h"

namespace feasst {

FEASST_MAPPER(Scattering,);

Scattering::Scattering(argtype * args) : Analyze(args) {
  num_frequency_ = integer("num_frequency", args, 100);
}
Scattering::Scattering(argtype args) : Scattering(&args) {
  feasst_check_all_used(args);
}

void Scattering::initialize(MonteCarlo * mc) {
  const Configuration& config = configuration(mc->system());
  ASSERT(config.dimension() == 3, "only implemented for 3d.");
  ASSERT(!config.domain().is_tilted(), "not implemented for tilted domains.");
  Position kvec;
  kvecs_.clear();
  site_ff_.clear();
  for (int kx = -num_frequency_; kx < num_frequency_; ++kx) {
    for (int ky = -num_frequency_; ky < num_frequency_; ++ky) {
      for (int kz = -num_frequency_; kz < num_frequency_; ++kz) {
        if (kx != 0 || ky != 0 || kz != 0) {
          kvec.set_vector({2*PI*kx/config.domain().side_length(0),
                           2*PI*ky/config.domain().side_length(1),
                           2*PI*kz/config.domain().side_length(2)});
          kvecs_.push_back(kvec);
          const double k = kvec.distance();
          std::vector<double> ff;
          for (int site_type = 0; site_type < config.num_site_types(); ++site_type) {
            const double radius = 0.5*config.model_params().select("sigma").value(site_type);
            const double volume = 4*PI*std::pow(radius, 3)/3.;
            ff.push_back(volume*3*(std::sin(k*radius)-k*radius*std::cos(k*radius))/std::pow(k*radius, 3));
          }
          site_ff_.push_back(ff);
        }
      }
    }
  }
  iq_.clear();
  iq_.resize(num_vectors());
}

void Scattering::update(const MonteCarlo& mc) {
  std::vector<double> fq(2*num_vectors());
  const Configuration& config = configuration(mc.system());
  const Select& selection = config.group_select(0);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        for (int k = 0; k < num_vectors(); ++k) {
          const double kr = site.position().dot_product(kvecs_[k]);
          const double ff = site_ff_[k][site.type()];
          fq[k] += ff*std::cos(kr);
          fq[k+num_vectors()] -= ff*std::sin(kr);
        }
      }
    }
  }
  for (int k = 0; k < num_vectors(); ++k) {
    const double fqr = fq[k];
    const double fqi = fq[k + num_vectors()];
    iq_[k].accumulate((fqr*fqr + fqi*fqi)/config.num_sites());
  }
  //for (int k = 0; k < 2*num_vectors(); ++k) {
  //  fq_[k].accumulate(fq[k]);
  //}
}

//std::vector<double> Scattering::iq_() const {
//  std::vector<double> iq(num_vectors());
//  for (int k = 0; k < num_vectors(); ++k) {
//    iq[k] = std::pow(fq_[k].average(), 2) + std::pow(fq_[k + num_vectors()].average(), 2);
//  }
//  return iq;
//}

std::string Scattering::write(const MonteCarlo& mc) {
  const Configuration& config = configuration(mc.system());
  //std::vector<double> iq = iq_();
  const int num_site_types = config.num_site_types();
  std::stringstream ss;
  ss << "#\"num_site_types\":" << num_site_types << "," << std::endl;
  ss << "q,i,";
  for (int site_type = 0; site_type < num_site_types; ++site_type) {
    ss << "p" << site_type << ",";
  }
  ss << std::endl;
  for (int k = 0; k < num_vectors(); ++k) {
    ss << kvecs_[k].distance() << "," << iq_[k].average() << ",";
    for (int site_type = 0; site_type < num_site_types; ++site_type) {
      ss << site_ff_[k][site_type] << ",";
    }
    ss << std::endl;
  }
  return ss.str();
}

void Scattering::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(6301, ostr);
  feasst_serialize(num_frequency_, ostr);
}

Scattering::Scattering(std::istream& istr)
  : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(6301 == version, "version mismatch:" << version);
  feasst_deserialize(&num_frequency_, istr);
}

}  // namespace feasst
