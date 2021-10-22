#include <cmath>  // exp
#include <fstream>
#include <algorithm>
#include <iomanip>  // setprecision
#include "utils/include/progress_report.h"
#include "utils/include/utils.h"  // is_equal
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "math/include/accumulator.h"
#include "flat_histogram/include/transition_matrix.h"

namespace feasst {

TransitionMatrix::TransitionMatrix(argtype args) : TransitionMatrix(&args) {
  check_all_used(args);
}
TransitionMatrix::TransitionMatrix(argtype * args) {
  class_name_ = "TransitionMatrix";
  min_visits_ = integer("min_visits", args, 100);
  min_sweeps_ = integer("min_sweeps", args);
  reset_sweeps_ = integer("reset_sweeps", args, -1);
  collection_ = TripleBandedCollectionMatrix(args);
}

void TransitionMatrix::update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool revert) {
  DEBUG("macro old/new " << macrostate_old << " " << macrostate_new);
  DEBUG("is_accepted " << is_accepted);
  const int bin = bin_(macrostate_old, macrostate_new, is_accepted);
  const int index = macrostate_new - macrostate_old + 1;
  DEBUG("bin " << bin << " index " << index);
  ASSERT(index >= 0 and index <= 2, "index(" << index << ") must be 0, 1 or 2");
  if (is_accepted && (macrostate_old != macrostate_new)) {
    if (revert) {
      --visits_[bin];
    } else {
      ++visits_[bin];
    }
  }
  double metropolis_prob = std::min(1., std::exp(ln_metropolis_prob));
  double reverse_prob = 1. - metropolis_prob;
  if (revert) {
    metropolis_prob *= -1.;
    reverse_prob *= -1.;
  }
  DEBUG("macrostate_old " << macrostate_old << " index " << index);
  DEBUG("metropolis_prob " << metropolis_prob);
  collection_.increment(macrostate_old, index, metropolis_prob);
  collection_.increment(macrostate_old, 1, reverse_prob);
  TRACE("colmat " << feasst_str(collection_.matrix()));
}

void TransitionMatrix::resize(const Histogram& histogram) {
  ln_prob_.resize(histogram.size());
  visits_.resize(histogram.size());
  collection_.resize(histogram.size());
  DEBUG("sizing bias " << histogram.size());
}

std::string TransitionMatrix::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << "num_sweeps," << num_sweeps_ << std::endl;
  TRACE("matrix," << feasst_str(collection_.matrix()));
  return ss.str();
}

std::string TransitionMatrix::write_per_bin_header() const {
  std::stringstream ss;
  ss << Bias::write_per_bin_header() << ",";
  ss << "visits,";
  ss << collection_.write_per_bin_header();
  return ss.str();
}

std::string TransitionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << ",";
  ss << MAX_PRECISION << visits_[bin] << ",";
  ss << collection_.write_per_bin(bin);
  return ss.str();
}

void TransitionMatrix::infrequent_update() {
  DEBUG("TransitionMatrix::infrequent_update_()");
  DEBUG("update the macrostate distribution");
  collection_.compute_ln_prob(&ln_prob_);

  DEBUG("update the number of sweeps");
  if (*std::min_element(visits_.begin(), visits_.end()) >= min_visits_) {
    ++num_sweeps_;
    std::fill(visits_.begin(), visits_.end(), 0);
    if (num_sweeps_ == reset_sweeps_) {
       increment_phase();
    }
  }

  DEBUG("check if complete");
  if (num_sweeps_ >= min_sweeps_) {
    set_complete_();
  }
  DEBUG("here");
}

void TransitionMatrix::set_ln_prob(
    const LnProbability& ln_prob) {
  ASSERT(ln_prob.size() == ln_prob_.size(), "size mismatch: " <<
    ln_prob.size() << " " << ln_prob_.size());
  ln_prob_ = ln_prob;
}

class MapTransitionMatrix {
 public:
  MapTransitionMatrix() {
    auto obj = MakeTransitionMatrix({{"min_sweeps", "0"}});
    obj->deserialize_map()["TransitionMatrix"] = obj;
  }
};

static MapTransitionMatrix mapper_ = MapTransitionMatrix();

std::shared_ptr<Bias> TransitionMatrix::create(std::istream& istr) const {
  return std::make_shared<TransitionMatrix>(istr);
}

TransitionMatrix::TransitionMatrix(std::istream& istr)
  : Bias(istr) {
  ASSERT(class_name_ == "TransitionMatrix", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(667 == version || version == 668, "mismatch version: " << version);
  feasst_deserialize_fstobj(&ln_prob_, istr);
  feasst_deserialize_fstobj(&collection_, istr);
  feasst_deserialize(&visits_, istr);
  feasst_deserialize(&min_visits_, istr);
  feasst_deserialize(&num_sweeps_, istr);
  feasst_deserialize(&min_sweeps_, istr);
  feasst_deserialize(&reset_sweeps_, istr);
  if (version == 667) {  // HWH backwards compatability
    int num_blocks, iter_block;
    bool is_block;
    std::vector<TransitionMatrix> blocks;
    feasst_deserialize(&num_blocks, istr);
    feasst_deserialize(&is_block, istr);
    feasst_deserialize(&iter_block, istr);
    feasst_deserialize_fstobj(&blocks, istr);
  }
}

void TransitionMatrix::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(668, ostr);
  feasst_serialize_fstobj(ln_prob_, ostr);
  feasst_serialize_fstobj(collection_, ostr);
  feasst_serialize(visits_, ostr);
  feasst_serialize(min_visits_, ostr);
  feasst_serialize(num_sweeps_, ostr);
  feasst_serialize(min_sweeps_, ostr);
  feasst_serialize(reset_sweeps_, ostr);
}

bool TransitionMatrix::is_equal(const TransitionMatrix& transition_matrix,
    const double tolerance) const {
  if (!collection_.is_equal(transition_matrix.collection_, tolerance)) {
    return false;
  }
  if (!ln_prob_.is_equal(transition_matrix.ln_prob_, tolerance)) {
    return false;
  }
  if (!feasst::is_equal(visits_, transition_matrix.visits_)) {
    INFO("visits not equal");
    return false;
  }
  if (min_visits_ != transition_matrix.min_visits_) return false;
  if (num_sweeps_ != transition_matrix.num_sweeps_) return false;
  if (min_sweeps_ != transition_matrix.min_sweeps_) return false;
  if (reset_sweeps_ != transition_matrix.reset_sweeps_) return false;
  return true;
}

void TransitionMatrix::set_num_iterations_to_complete(const int sweeps) {
  min_sweeps_ = sweeps;
  if (num_sweeps_ < min_sweeps_) set_incomplete_();
}

TransitionMatrix::TransitionMatrix(const Bias& bias) {
  std::stringstream ss;
  bias.serialize(ss);
  *this = TransitionMatrix(ss);
}

}  // namespace feasst
