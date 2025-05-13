
#include <algorithm>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "math/include/histogram.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/transition_matrix_guess.h"

namespace feasst {

FEASST_MAPPER(TransitionMatrixGuess, argtype({{"min_sweeps", "1"}}));

TransitionMatrixGuess::TransitionMatrixGuess(argtype * args) {
  class_name_ = "TransitionMatrixGuess";
  min_collect_sweeps_ = integer("min_collect_sweeps", args, -1);
  const std::string ln_prob_file = str("ln_prob_file", args, "");
  if (!ln_prob_file.empty()) {
    ln_prob_guess_ = std::make_unique<LnProbability>();
    ln_prob_guess_->read(ln_prob_file);
  }
  transition_matrix_ = std::make_unique<TransitionMatrix>(args);
}
TransitionMatrixGuess::TransitionMatrixGuess(argtype args) : TransitionMatrixGuess(&args) {
  feasst_check_all_used(args);
}
TransitionMatrixGuess::~TransitionMatrixGuess() {}

bool TransitionMatrixGuess::is_tm_bias_(const Macrostate& macro) const {
  if (ln_prob_guess_ &&
      transition_matrix_->num_cycles(-1, macro) < min_collect_sweeps_) {
    return false;
  }
  return true;
}

void TransitionMatrixGuess::update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const Macrostate& macro) {
  transition_matrix_->update(macrostate_old, macrostate_new,
    ln_metropolis_prob, is_accepted, is_endpoint, macro);
  is_tm_bias_at_update_ = is_tm_bias_(macro);
  if (is_tm_bias_at_update_) {
    if (production_ == 0) {
      production_ = 1;
      increment_phase();
    }
    if (transition_matrix_->is_complete()) {
      set_complete_();
    }
  }
}

const LnProbability& TransitionMatrixGuess::ln_prob() const {
  if (is_tm_bias_at_update_) {
    return transition_matrix_->ln_prob();
  } else {
    ASSERT(ln_prob_guess_->size() == transition_matrix_->ln_prob().size(),
      "The LnProbability guess size: " << ln_prob_guess_->size() << "is not "
      "equal to the TransitionMatrix size: "
      << transition_matrix_->ln_prob().size());
    return *ln_prob_guess_;
  }
}

void TransitionMatrixGuess::resize(const Histogram& histogram) {
  ln_prob_guess_->resize(histogram.size());
  transition_matrix_->resize(histogram);
}

void TransitionMatrixGuess::infrequent_update(const Macrostate& macro) {
  transition_matrix_->infrequent_update(macro);
}

std::string TransitionMatrixGuess::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << transition_matrix_->write();
  return ss.str();
}

std::string TransitionMatrixGuess::write_per_bin_header(const std::string& append) const {
  std::stringstream ss;
  ss << transition_matrix_->write_per_bin_header("");
  return ss.str();
}

std::string TransitionMatrixGuess::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << transition_matrix_->write_per_bin(bin);
  return ss.str();
}

void TransitionMatrixGuess::set_ln_prob(
    const LnProbability& ln_prob) {
  FATAL("not implemented. Try WangLandau or TransitionMatrix for reweighting.");
}

std::shared_ptr<Bias> TransitionMatrixGuess::create(std::istream& istr) const {
  return std::make_shared<TransitionMatrixGuess>(istr);
}

int TransitionMatrixGuess::num_cycles(const int state, const Macrostate& macro) const {
  return transition_matrix_->num_cycles(state, macro);
}

bool TransitionMatrixGuess::is_adjust_allowed(const Macrostate& macro) const {
  FATAL("not implemented");
}

TransitionMatrixGuess::TransitionMatrixGuess(std::istream& istr) : Bias(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6327, "mismatch version: " << version);
  feasst_deserialize(&min_collect_sweeps_, istr);
  feasst_deserialize(&production_, istr);
  feasst_deserialize(ln_prob_guess_, istr);
  feasst_deserialize(transition_matrix_, istr);
}

void TransitionMatrixGuess::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(6327, ostr);
  feasst_serialize(min_collect_sweeps_, ostr);
  feasst_serialize(production_, ostr);
  feasst_serialize(ln_prob_guess_, ostr);
  feasst_serialize(transition_matrix_, ostr);
}

int TransitionMatrixGuess::cycles_to_complete() const {
  return transition_matrix_->cycles_to_complete(); }

void TransitionMatrixGuess::set_cycles_to_complete(const int sweeps) {
  transition_matrix_->set_cycles_to_complete(sweeps); }

const TransitionMatrix& TransitionMatrixGuess::transition_matrix() const {
  return const_cast<TransitionMatrix&>(*transition_matrix_); }

void TransitionMatrixGuess::set_cm(const int macro, const Bias& bias) {
  transition_matrix_->set_cm(macro, bias); }

const CollectionMatrix& TransitionMatrixGuess::cm() const {
  return transition_matrix().cm(); }

const int TransitionMatrixGuess::visits(const int macro, const int index) const {
  return transition_matrix().visits(macro, index); }

}  // namespace feasst
