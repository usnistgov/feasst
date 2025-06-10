#include <cmath>
#include "utils/include/arguments.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle.h"
#include "configuration/include/site.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/monte_carlo.h"
#include "fftw/include/scattering_fftw.h"

namespace feasst {

FEASST_MAPPER(ScatteringFFTW,);

ScatteringFFTW::ScatteringFFTW(argtype * args) : Analyze(args) {
  bin_spacing_ = dble("bin_spacing", args, 0.1);
  bins_per_side_ = integer("bins_per_side", args, -1);
  delta_rho_ = dble("delta_rho", args, 1);
}
ScatteringFFTW::ScatteringFFTW(argtype args) : ScatteringFFTW(&args) {
  feasst_check_all_used(args);
}

void ScatteringFFTW::resize_fftw_variables_() {
  in_  = reinterpret_cast<double*>(fftw_malloc(sizeof(double) * feasst::product(num_bin_)));
  DEBUG("size of in: " << feasst::product(num_bin_));
  out_ = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * num_q_()));
  plan_ = fftw_plan_dft_r2c_3d(num_bin_[0], num_bin_[1], num_bin_[2], in_, out_, FFTW_MEASURE);
  fftw_initialized_ = true;
}

void ScatteringFFTW::initialize(MonteCarlo * mc) {
  const Configuration& config = configuration(mc->system());
  ASSERT(config.dimension() == 3, "only implemented for 3d.");
  ASSERT(!config.domain().is_tilted(), "not implemented for tilted domains.");

  // fftw
  updates_ = 0;
  num_bin_.resize(config.dimension());
  bin_dist_.resize(config.dimension());
  lower_bound_.resize(config.dimension());
  DEBUG("bin_spacing " << bin_spacing_);
  for (int dim = 0; dim < config.dimension(); ++dim) {
    const double side_length = config.domain().side_length(dim);
    DEBUG("side_length " << side_length);
    if (bins_per_side_ > 0) {
      ASSERT(config.domain().is_cubic(), "bins_per_side assumes cubic Domain.");
      num_bin_[dim] = bins_per_side_;
    } else {
      num_bin_[dim] = static_cast<int>(side_length/bin_spacing_);
      if (num_bin_[dim] % 2 != 0) ++num_bin_[dim];
    }
    bin_dist_[dim] = side_length/static_cast<double>(num_bin_[dim]);
    lower_bound_[dim] = -side_length/2.;
  }
  DEBUG("num_bin " << feasst_str(num_bin_));
  DEBUG("bin_dist " << feasst_str(bin_dist_));
  sqsq_.resize(num_q_());
  fksq_.resize(num_q_());
  resize_fftw_variables_();
}

void ScatteringFFTW::fill_grid_(const Configuration& config, const bool centers) {
  const int nx = num_bin_[0];
  const int ny = num_bin_[1];
  const int nz = num_bin_[2];
  // zero grid
  for (int ix = 0; ix < nx; ++ix) {
  for (int iy = 0; iy < ny; ++iy) {
  for (int iz = 0; iz < nz; ++iz) {
    in_[iz + nz*(iy + ny*ix)] = 0;
  }}}

  std::vector<int> center(config.dimension());
  const Select& selection = config.group_select(0);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        for (int dim = 0; dim < config.dimension(); ++dim) {
          const int num_bin = num_bin_[dim];
          int * cen = &center[dim];
          *cen = static_cast<int>((site.position().coord(dim) - lower_bound_[dim])/bin_dist_[dim]);
          // wrap positions in grid
          if (*cen >= num_bin) {
            *cen -= num_bin;
          } else if (*cen < 0) {
            *cen += num_bin;
          }
          ASSERT((*cen < num_bin) && (*cen >= 0),
            "Not wrapped: " << *cen << " num_bin:" << num_bin);
        }
        DEBUG("pos " << site.position().str());
        DEBUG("center " << feasst_str(center));
        if (centers) {
          const int bin = center[2] + nz*(center[1] + ny*center[0]);
          DEBUG("bin " << bin);
          ASSERT(bin < feasst::product(num_bin_), "Err. bin:" << bin <<
            " num:" << feasst::product(num_bin_));
          in_[bin] = std::pow(-1, center[0] + center[1] + center[2]);
        } else {
          FATAL("not implemented");
//          // implement a recursive search for which bins overlap with the site
//          if (site.is_anisotropic()) {
//            FATAL("not implemented");
//          } else {
//            const double sigma = config.model_params().select("sigma").value(site.type());
//            fill_grid_site_(sigma, site.position(), center);
//          }
        }
      }
    }
  }
}

void ScatteringFFTW::update(const MonteCarlo& mc) {
  const Configuration& config = configuration(mc.system());

  // sq
  fill_grid_(config, true);
  fftw_execute(plan_);
  for (int iq = 0; iq < num_q_(); ++iq) {
    sqsq_[iq] += (std::pow(out_[iq][0], 2) + std::pow(out_[iq][1], 2))/
                 static_cast<double>(config.num_particles());
  }

//  // iq
//  fill_grid_(config, false);
//  fftw_execute(plan_);
//  for (int iq = 0; iq < num_q_(); ++iq) {
//    fksq_[iq] += std::pow(out_[iq][0], 2) + std::pow(out_[iq][1], 2);
//  }

  ++updates_;
}

std::string ScatteringFFTW::write(const MonteCarlo& mc) {
  const Configuration& config = configuration(mc.system());
  // integerate 3D->1D
  std::vector<double> sq1d, fk1d, count;
  const int nx = num_bin_[0];
  const int ny = num_bin_[1];
  const int nz = num_bin_[2]/2 + 1;
  const double cx = nx/2.;
  const double cy = ny/2.;
  const double cz = num_bin_[2]/2.;
  for (int ix = 0; ix < nx; ++ix) {
  for (int iy = 0; iy < ny; ++iy) {
  for (int iz = 0; iz < nz; ++iz) {
    const double dx = static_cast<double>(ix) - cx;
    const double dy = static_cast<double>(iy) - cy;
    const double dz = static_cast<double>(iz) - cz;
    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double round = r/delta_rho_;
    const int lower = static_cast<int>(round);
    const int upper = lower + 1;
    if (upper >= static_cast<int>(sq1d.size())) {
      sq1d.resize(upper + 1);
      fk1d.resize(upper + 1);
      count.resize(upper + 1);
    }
    const int index3d = iz + nz*(iy + ny*ix);
    const double ufrac = round - static_cast<double>(lower);
    sq1d[upper] += ufrac*sqsq_[index3d];
    fk1d[upper] += ufrac*fksq_[index3d];
    count[upper] += ufrac;
    sq1d[lower] += (1. - ufrac)*sqsq_[index3d];
    fk1d[lower] += (1. - ufrac)*fksq_[index3d];
    count[lower] += 1. - ufrac;
  }}}

  // normalize by counts and updates
  for (int iq = 0; iq < static_cast<int>(sq1d.size()); ++iq) {
    sq1d[iq] /= count[iq] * updates_;
    fk1d[iq] /= count[iq] * updates_;
  }

  // compute the minimum index, which is the smallest q value of 2pi/L
  int min_index = -1;
  { int index = 0;
    while(min_index < 0) {
      if (delta_rho_*static_cast<double>(index) >= 1. - 100*NEAR_ZERO) {
        min_index = index;
      }
      ++index;
    }
  }

  // compute invariants to normalize intensity
//  double vp = 1;
//  const double num_sites = static_cast<double>(config.num_sites());
  const double volume = config.domain().volume();
//  const double phi = num_sites*vp/volume;
//  const double invariant = 2*PI*PI*phi*(1-phi);
  const double fk1dmin = *std::min_element(fk1d.begin()+min_index, fk1d.end());
  std::vector<double> integrand;
  for (int iq = min_index; iq < static_cast<int>(fk1d.size()); ++iq) {
    const double q = delta_rho_*static_cast<double>(iq)*2.*PI/std::pow(volume, 1./3.);
    integrand.push_back(q*q*(fk1d[iq] - fk1dmin));
  }

  std::stringstream ss;
//  ss << "#\"num_site_types\":" << num_site_types << "," << std::endl;
  ss << "q,i,count,index" << std::endl;
  //ss << "q,i,in,s,count,index" << std::endl;
  for (int iq = min_index; iq < static_cast<int>(sq1d.size()); ++iq) {
    const double q = delta_rho_*static_cast<double>(iq)*2.*PI/std::pow(volume, 1./3.);
    ss << q << ","
//       << fk1d[iq]/num_sites << ",";
//    fk1d[iq] *= invariant * 1e8;
//    ss << fk1d[iq] << ","
       << sq1d[iq] << ","
       << count[iq] << ","
       << iq << std::endl;
  }
  return ss.str();
}

ScatteringFFTW::~ScatteringFFTW() {
  if (fftw_initialized_) {
    fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
  }
}

void ScatteringFFTW::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(3689, ostr);
  feasst_serialize(fftw_initialized_, ostr);
  feasst_serialize(bin_spacing_, ostr);
  feasst_serialize(bins_per_side_, ostr);
  feasst_serialize(delta_rho_, ostr);
  feasst_serialize(updates_, ostr);
  feasst_serialize(num_bin_, ostr);
  feasst_serialize(bin_dist_, ostr);
  feasst_serialize(lower_bound_, ostr);
  feasst_serialize(fksq_, ostr);
  feasst_serialize(sqsq_, ostr);
}

ScatteringFFTW::ScatteringFFTW(std::istream& istr)
  : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(3689 == version, "version mismatch:" << version);
  feasst_deserialize(&fftw_initialized_, istr);
  feasst_deserialize(&bin_spacing_, istr);
  feasst_deserialize(&bins_per_side_, istr);
  feasst_deserialize(&delta_rho_, istr);
  feasst_deserialize(&updates_, istr);
  feasst_deserialize(&num_bin_, istr);
  feasst_deserialize(&bin_dist_, istr);
  feasst_deserialize(&lower_bound_, istr);
  feasst_deserialize(&fksq_, istr);
  feasst_deserialize(&sqsq_, istr);
  if (fftw_initialized_) {
    resize_fftw_variables_();
  }
}

}  // namespace feasst
