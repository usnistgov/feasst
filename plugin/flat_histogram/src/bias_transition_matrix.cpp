
#include <algorithm>
#include "flat_histogram/include/bias_transition_matrix.h"
#include "core/include/utils_io.h"
#include "core/include/utils_math.h"
#include "core/include/debug.h"
#include "core/include/accumulator.h"

namespace feasst {

BiasTransitionMatrix::BiasTransitionMatrix(const argtype &args) {
  args_.init(args);
  min_visits_ = args_.key("min_visits").dflt("100").integer();
  min_sweeps_ = args_.key("min_sweeps").integer();
  num_steps_to_update_ =
    args_.key("num_steps_to_update").dflt("1000000").integer();
}

void BiasTransitionMatrix::update(const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted) {
  DEBUG("macro old/new " << macrostate_old << " " << macrostate_new);
  const int bin = bin_(macrostate_old, macrostate_new, is_accepted);
  const int index = macrostate_new - macrostate_old + 1;
  DEBUG("bin " << bin << " index " << index);
  ASSERT(index >= 0 && index <= 2, "index(" << index << ") must be 0, 1 or 2");
  if (is_accepted) {
    ++visits_[bin];
  }
  const double metropolis_prob = std::min(1., exp(ln_metropolis_prob));
  DEBUG("macrostate_old " << macrostate_old << " index " << index);
  collection_.increment(macrostate_old, index, metropolis_prob);
  collection_.increment(macrostate_old, 1, 1. - metropolis_prob);
  update_blocks_(macrostate_old, macrostate_new,
                 ln_metropolis_prob, is_accepted);
  ++num_steps_since_update_;
  DEBUG("num_steps_since_update_ " << num_steps_since_update_ << " num_steps_to_update_ " << num_steps_to_update_);
  if (num_steps_since_update_ >= num_steps_to_update_) {
    num_steps_since_update_ = 0;
    DEBUG("updating");
    infrequent_update_();
  }
}

void BiasTransitionMatrix::update_blocks_(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted) {
  // If the object that is updating is a block, don't update blocks or you'll
  // have an infinite recursion.
  if (is_block_) {
    return;
  }

  // If the blocks haven't been initialized, create new blocks.
  if (blocks_.size() == 0) {
    for (int index = 0; index < 3; ++index) {
      auto block = std::make_shared<BiasTransitionMatrix>(*this);
      block->is_block_ = true;
      block->blocks_.clear();
      blocks_.push_back(block);
    }
  }

  // Update one of the randomly chosen blocks.
  random_.element(blocks_)->update(macrostate_old,
                                   macrostate_new,
                                   ln_metropolis_prob,
                                   is_accepted);
}

void BiasTransitionMatrix::resize(const Histogram& histogram) {
  ln_macro_prob_.resize(histogram.size());
  visits_.resize(histogram.size());
  collection_.resize(histogram.size());
  DEBUG("sizing bias " << histogram.size());
}

std::string BiasTransitionMatrix::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << "num_sweeps " << num_sweeps_ << std::endl;
  DEBUG("matrix " << feasst_str(collection_.matrix()));
  return ss.str();
}

std::string BiasTransitionMatrix::write_per_bin_header() const {
  std::stringstream ss;
  ss << Bias::write_per_bin_header() << " ";
  for (const std::shared_ptr<BiasTransitionMatrix> block : blocks_) {
    ss << "lnpi_partial ";
   // block->write_per_bin_header();
  }
  if (!is_block_) {
    ss << "lnpi_stdev ";
  }
  ss << "c0 c1 c2";
  return ss.str();
}

std::string BiasTransitionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << " ";

  // compute and print the standard deviation of the average of the blocks
  Accumulator acc_block;
  for (const std::shared_ptr<BiasTransitionMatrix> block : blocks_) {
    ss << block->ln_macro_prob().value(bin) << " ";
   // write_per_bin(bin) << " ";
    DEBUG(block->ln_macro_prob().value(bin));
    acc_block.accumulate(block->ln_macro_prob().value(bin));
  }
  if (!is_block_) {
    ss << acc_block.stdev_of_av() << " ";
//    const std::vector<double>& cols =  collection_.matrix()[bin];
//    ss << "sz " << cols.size() << " " ;
//    for (auto element : cols) {
//    //for (auto element : collection_.matrix()[bin]) {
//      ss << element << " ";
//    }
//    ss << "size " << collection_.matrix()[bin].size() << " ";
    ss << collection_.matrix()[bin][0] << " "
       << collection_.matrix()[bin][1] << " "
       << collection_.matrix()[bin][2];
//    ss << " " << cols[0] << " "
//      << cols[1] << " "
//      << cols[2];
  }
  return ss.str();
}

void BiasTransitionMatrix::infrequent_update_() {
  DEBUG("BiasTransitionMatrix::infrequent_update_() " << is_block_);
  // update the macrostate distribution
  collection_.compute_ln_prob(&ln_macro_prob_);
  for (const std::shared_ptr<BiasTransitionMatrix> block : blocks_) {
    block->infrequent_update_();
  }

  // update the number of sweeps
  if (*std::min_element(visits_.begin(), visits_.end()) >= min_visits_) {
    ++num_sweeps_;
    std::fill(visits_.begin(), visits_.end(), 0);
  }

  if (num_sweeps_ >= min_sweeps_) {
    set_complete_();
  }
}

}  // namespace feasst
