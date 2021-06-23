#include <cmath>  // exp
#include <fstream>
#include <algorithm>
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
  num_blocks_ = integer("num_blocks", args, 30);
  reset_sweeps_ = integer("reset_sweeps", args, -1);
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
  update_blocks_(macrostate_old, macrostate_new,
                 ln_metropolis_prob, is_accepted, revert);
  if (!is_block_ && !dump_file_.empty()) {
    std::ofstream file;
    file.open(dump_file_, std::ofstream::out | std::ofstream::app);
    file << macrostate_old << " " << macrostate_new << " "
         //<< ln_metropolis_prob << std::endl;
         << MAX_PRECISION << ln_metropolis_prob << std::endl;
  }
}

void TransitionMatrix::update_blocks_(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool revert) {
  // If the object that is updating is a block, don't update blocks or you'll
  // have an infinite recursion.
  if (is_block_) return;

  // If the blocks haven't been initialized, create new blocks.
  if (blocks_.size() == 0) {
    for (int index = 0; index < num_blocks_; ++index) {
      TransitionMatrix block(*this);
      block.is_block_ = true;
      block.blocks_.clear();
      blocks_.push_back(block);
    }
  }

  // Update one of the blocks.
  if (revert) {
    if (iter_block_ == 0) {
      iter_block_ = static_cast<int>(blocks_.size()) - 1;
    } else {
      --iter_block_;
    }
  } else {
    ++iter_block_;
    if (iter_block_ == static_cast<int>(blocks_.size())) {
      iter_block_ = 0;
    }
  }
  if (blocks_.size() > 0) {
    blocks_[iter_block_].update_or_revert(macrostate_old,
                                          macrostate_new,
                                          ln_metropolis_prob,
                                          is_accepted,
                                          revert);
  }
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
//  for (int index = 0; index < static_cast<int>(blocks_.size()); ++index) {
//  // for (const TransitionMatrix& block : blocks_) {
//    ss << "ln_prob_" << index << ",";
//   // block->write_per_bin_header();
//  }
  if (blocks_.size() > 2) ss << "ln_prob_stdev,";
  ss << "visits,c0,c1,c2,";
  return ss.str();
}

std::string TransitionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << ",";

  // compute and print the standard deviation of the average of the blocks
  Accumulator acc_block;
  for (const TransitionMatrix& block : blocks_) {
//    ss << MAX_PRECISION << block.ln_prob().value(bin) << ",";
   // write_per_bin(bin) << ",";
    TRACE(block.ln_prob().value(bin));
    acc_block.accumulate(block.ln_prob().value(bin));
  }
  if (blocks_.size() > 2) {
    ss << MAX_PRECISION << acc_block.stdev_of_av() << ",";
//    const std::vector<double>& cols =  collection_.matrix()[bin];
//    ss << "sz " << cols.size() << " " ;
//    for (auto element : cols) {
//    //for (auto element : collection_.matrix()[bin]) {
//      ss << element << " ";
//    }
//    ss << "size " << collection_.matrix()[bin].size() << " ";
  }
  ss << MAX_PRECISION
     << visits_[bin] << ","
     << collection_.matrix()[bin][0] << ","
     << collection_.matrix()[bin][1] << ","
     << collection_.matrix()[bin][2] << ",";
//    ss << " " << cols[0] << " "
//      << cols[1] << " "
//      << cols[2];
  return ss.str();
}

void TransitionMatrix::infrequent_update() {
  DEBUG("TransitionMatrix::infrequent_update_() " << is_block_);
  DEBUG("update the macrostate distribution");
  collection_.compute_ln_prob(&ln_prob_);

  DEBUG("update the blocks");
  for (TransitionMatrix& block : blocks_) {
    block.infrequent_update();
  }

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
  ASSERT(667 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&ln_prob_, istr);
  feasst_deserialize_fstobj(&collection_, istr);
  feasst_deserialize(&visits_, istr);
  feasst_deserialize(&min_visits_, istr);
  feasst_deserialize(&num_sweeps_, istr);
  feasst_deserialize(&min_sweeps_, istr);
  feasst_deserialize(&reset_sweeps_, istr);
  feasst_deserialize(&num_blocks_, istr);
  feasst_deserialize(&is_block_, istr);
  feasst_deserialize(&iter_block_, istr);
  feasst_deserialize_fstobj(&blocks_, istr);
}

void TransitionMatrix::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(667, ostr);
  feasst_serialize_fstobj(ln_prob_, ostr);
  feasst_serialize_fstobj(collection_, ostr);
  feasst_serialize(visits_, ostr);
  feasst_serialize(min_visits_, ostr);
  feasst_serialize(num_sweeps_, ostr);
  feasst_serialize(min_sweeps_, ostr);
  feasst_serialize(reset_sweeps_, ostr);
  feasst_serialize(num_blocks_, ostr);
  feasst_serialize(is_block_, ostr);
  feasst_serialize(iter_block_, ostr);
  feasst_serialize_fstobj(blocks_, ostr);
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
  if (num_blocks_ != transition_matrix.num_blocks_) return false;
  return true;
}

void TransitionMatrix::set_num_iterations(const int sweeps) {
  min_sweeps_ = sweeps;
  if (num_sweeps_ < min_sweeps_) set_incomplete_();
}

TransitionMatrix::TransitionMatrix(const Bias& bias) {
  std::stringstream ss;
  bias.serialize(ss);
  *this = TransitionMatrix(ss);
}

void TransitionMatrix::read_dump_file(const std::string dump_file) {
  std::ifstream file(dump_file);
  int macro_old, macro_new;
  double ln_prob;
  ProgressReport prog({{"num", feasst::str(1e8)}});
  while (!file.eof()) {
    std::string line;
    file >> macro_old >> macro_new >> ln_prob;
    if (!file.eof()) {
      DEBUG("old " << macro_old << " new " << macro_new << " lnp " << ln_prob);
      update_or_revert(macro_old, macro_new, ln_prob, false, false);
      prog.check();
    }
  }
}

}  // namespace feasst
