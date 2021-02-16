#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/ensemble.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/clones.h"

namespace feasst {

void Ensemble::init_(const Histogram& macrostates,
                     const LnProbability& ln_prob) {
  ln_prob_original_ = ln_prob;
  macrostates_ = macrostates;
  ASSERT(ln_prob_original().size() == macrostates_.size(),
    "ln_prob(" << ln_prob_original().size() << ") or macrostate (" <<
    macrostates_.size() << " ) are not of same size.");
  ln_prob_ = ln_prob_original();
}

Ensemble::Ensemble(const FlatHistogram& flat_hist)
  : Ensemble(flat_hist.macrostate().histogram(),
             flat_hist.ln_prob()) {}

Ensemble::Ensemble(const Clones& clones) {
  Histogram macrostates;
  LnProbability ln_prob = clones.ln_prob(&macrostates);
  init_(macrostates, ln_prob);
}

void Ensemble::phase_boundary(const int phase, int * min, int * max) const {
  std::vector<int> mins = ln_prob().minima();
  const double num_min = static_cast<int>(mins.size());
  if (num_min == 0) {
    *min = 0;
    *max = ln_prob().size() - 1;
  } else if (num_min == 1) {
    if (phase == 0) {
      *min = 0;
      *max = mins[0];
    } else if (phase == 1) {
      *min = mins[0];
      *max = ln_prob().size() - 1;
    } else {
      ERROR("unrecognized phase: " << phase);
    }
  } else {
    FATAL("multiple minima: " << num_min << " not implemented");
  }
}

bool Ensemble::is_phase_boundary() const {
  int min, max;
  phase_boundary(0, &min, &max);
  if (min == 0 && max == ln_prob().size() - 1) {
    return false;
  }
  return true;
}

const LnProbability& Ensemble::reweight(const double delta_conjugate) {
  delta_conjugate_ = delta_conjugate;
  ln_prob_ = ln_prob_original();
  for (int macro = 0; macro < ln_prob_.size(); ++macro) {
    ln_prob_.add(macro, macrostates().center_of_bin(macro)
                 *delta_conjugate);
  }
  ln_prob_.normalize();
  return ln_prob_;
}

double Ensemble::average(const std::vector<double>& macrostate_averages,
     const int phase) const {
  ASSERT(ln_prob().size() == static_cast<int>(macrostate_averages.size()),
    "size mismatch: ln_prob:" << ln_prob().size() <<
    " macro:" << macrostate_averages.size());
  int min, max;
  phase_boundary(phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostate_averages[bin]*std::exp(ln_prob().value(bin));
  }
  return average/ln_prob().sum_probability(min, max);
}

double Ensemble::average_macrostate(const int phase) const {
  int min, max;
  phase_boundary(phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostates().center_of_bin(bin)*std::exp(ln_prob().value(bin));
  }
  return average/ln_prob().sum_probability(min, max);
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(const Histogram& macrostates,
  const LnProbability& ln_prob,
  const double conjugate) : Ensemble(macrostates, ln_prob) {
  original_conjugate_ = conjugate;
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(
    const FlatHistogram& flat_histogram,
    const double beta_mu) : Ensemble(flat_histogram) {
  original_conjugate_ = beta_mu;
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(
    const Clones& clones) : Ensemble(clones) {
  original_conjugate_ = clones.clone(0).system().thermo_params().beta_mu();
}

double GrandCanonicalEnsemble::betaPV(const int phase) const {
  int min, max;
  phase_boundary(phase, &min, &max);
  return -ln_prob().value(0) + std::log(ln_prob().sum_probability(min, max));
}

void ExtrapolateBetaGCE::extrapolateBetaGCE_(
    const std::vector<std::vector<double> >& energy_moments,
    argtype args) {
  const int order = integer("order", &args, 2);
  ASSERT(order == 2, "Only second order is currently implemented");
  ASSERT(static_cast<int>(energy_moments.size()) == order,
    "order: " << order << " doesn't match energy: " << energy_moments.size());
  ASSERT(static_cast<int>(energy_moments[0].size()) == ln_prob_original_.size(),
    "size of energy:" << energy_moments[0].size() << " doesn't match ln_prob: " <<
    ln_prob_original_.size());
  ASSERT(!is_phase_boundary(), "assumes no phase boundary when averaging.");
  // add option to ignore phases, phase=-1
  ASSERT(std::abs(macrostates().center_of_bin(0)) < NEAR_ZERO,
    "assumes first marcostate is 0 particles for N and <NU>");
  const double beta_new = dble("beta_new", &args);
  const double beta_original = dble("beta_original", &args);
  check_all_used(args);
  const double dbeta = beta_new - beta_original;
  DEBUG("dbeta " << dbeta);
  std::vector<double> nu(ln_prob_original_.size());
  for (int bin = 0; bin < ln_prob_original_.size(); ++bin) {
    const double n = macrostates().center_of_bin(bin);
    nu[bin] = n*energy_moments[0][bin];
  }
  const double gc_u  = average(energy_moments[0]);
  const double gc_u2 = average(energy_moments[1]);
  const double gc_n = average_macrostate();
  const double gc_nu = average(nu);
  const double mu = original_conjugate()/beta_original;
  DEBUG("mu " << mu);
  DEBUG("gc_u " << gc_u);
  DEBUG("gc_u2 " << gc_u2);
  DEBUG("gc_n " << gc_n);
  DEBUG("gc_nu " << gc_nu);
  energy_ = energy_moments[0];
  for (int state = 0; state < ln_prob_original_.size(); ++state) {
    const double n = macrostates().center_of_bin(state);
    DEBUG("n " << n);
    const double u  = energy_moments[0][state];
    DEBUG("u " << u);
    const double u2 = energy_moments[1][state];
    DEBUG("u2 " << u2);
    const double dlnpdb = -u + gc_u + mu*n;
    DEBUG("dlnpdb " << dlnpdb);
    const double dudb = -u2 + u*u;
    DEBUG("dudb " << dudb);
    const double d2lnpdb2 = -dudb - gc_u2 + gc_u*gc_u - mu*(gc_nu - gc_n*gc_u);
    DEBUG("d2lnpdb2 " << d2lnpdb2);
    energy_[state] += dudb*dbeta;
    ln_prob_original_.add(state, dlnpdb*dbeta + d2lnpdb2*dbeta*dbeta/2.);
  }
}

ExtrapolateBetaGCE::ExtrapolateBetaGCE(
    const MonteCarlo& mc,
    const FlatHistogram& flat_hist,
    argtype args)
  : GrandCanonicalEnsemble(flat_hist, mc.system().thermo_params().beta_mu()) {
  const Analyze& an = SeekAnalyze().reference("Energy", mc);
  const int num_moments = static_cast<int>(an.accumulator().moments().size());
  DEBUG("num_moments " << num_moments);
  std::vector<std::vector<double> > energy_moments(num_moments);
  for (int mom = 0; mom < num_moments; ++mom) {
    energy_moments[mom] =
      SeekAnalyze().multistate_data("Energy", mc, AccumulatorMoment(mom));
  }
  extrapolateBetaGCE_(energy_moments, args);
}

ExtrapolateBetaGCE::ExtrapolateBetaGCE(const Clones& clones,
    argtype args) : GrandCanonicalEnsemble(clones) {
  const Analyze& an = SeekAnalyze().reference("Energy", clones.clone(0));
  const int num_moments = static_cast<int>(an.accumulator().moments().size());
  DEBUG("num_moments " << num_moments);
  std::vector<std::vector<double> > energy_moments(num_moments);
  for (int mom = 0; mom < num_moments; ++mom) {
    clones.stitch(&energy_moments[mom], "Energy", AccumulatorMoment(mom));
  }
  extrapolateBetaGCE_(energy_moments, args);
}

}  // namespace feasst
