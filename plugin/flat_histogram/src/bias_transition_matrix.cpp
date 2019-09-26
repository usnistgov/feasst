
#include <algorithm>
#include <math.h>
#include "flat_histogram/include/bias_transition_matrix.h"
#include "utils/include/utils_io.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"
#include "math/include/accumulator.h"

namespace feasst {

BiasTransitionMatrix::BiasTransitionMatrix(const argtype &args) {
  class_name_ = "BiasTransitionMatrix";
  Arguments args_(args);
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
  ASSERT(index >= 0 and index <= 2, "index(" << index << ") must be 0, 1 or 2");
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

void BiasTransitionMatrix::revert(const int macrostate_new,
    const int macrostate_old) {
  ERROR("not implemented: what if infrequent update needs to be reverted?");
}

void BiasTransitionMatrix::update_blocks_(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted) {
  // If the object that is updating is a block, don't update blocks or you'll
  // have an infinite recursion.
  if (is_block_) return;

  // If the blocks haven't been initialized, create new blocks.
  if (blocks_.size() == 0) {
    for (int index = 0; index < 3; ++index) {
      BiasTransitionMatrix block(*this);
      block.is_block_ = true;
      block.blocks_.clear();
      blocks_.push_back(block);
    }
  }

  // Update one of the blocks.
  ++iter_block_;
  if (iter_block_ == static_cast<int>(blocks_.size())) {
    iter_block_ = 0;
  }
  blocks_[iter_block_].update(macrostate_old,
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
  for (int index = 0; index < static_cast<int>(blocks_.size()); ++index) {
    ss << "lnpi_partial ";
   // block->write_per_bin_header();
  }
  if (!is_block_) ss << "lnpi_stdev ";
  ss << "c0 c1 c2";
  return ss.str();
}

std::string BiasTransitionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << " ";

  // compute and print the standard deviation of the average of the blocks
  Accumulator acc_block;
  for (const BiasTransitionMatrix& block : blocks_) {
    ss << block.ln_macro_prob().value(bin) << " ";
   // write_per_bin(bin) << " ";
    DEBUG(block.ln_macro_prob().value(bin));
    acc_block.accumulate(block.ln_macro_prob().value(bin));
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
  for (BiasTransitionMatrix& block : blocks_) {
    block.infrequent_update_();
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

class MapBiasTransitionMatrix {
 public:
  MapBiasTransitionMatrix() {
    auto obj = MakeBiasTransitionMatrix({{"min_sweeps", "0"}});
    obj->deserialize_map()["BiasTransitionMatrix"] = obj;
  }
};

static MapBiasTransitionMatrix mapper_ = MapBiasTransitionMatrix();

std::shared_ptr<Bias> BiasTransitionMatrix::create(std::istream& istr) const {
  return std::make_shared<BiasTransitionMatrix>(istr);
}

BiasTransitionMatrix::BiasTransitionMatrix(std::istream& istr)
  : Bias(istr) {
  ASSERT(class_name_ == "BiasTransitionMatrix", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(667 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&ln_macro_prob_, istr);
  feasst_deserialize_fstobj(&collection_, istr);
  feasst_deserialize(&visits_, istr);
  feasst_deserialize(&min_visits_, istr);
  feasst_deserialize(&num_sweeps_, istr);
  feasst_deserialize(&min_sweeps_, istr);
  feasst_deserialize(&num_steps_to_update_, istr);
  feasst_deserialize(&num_steps_since_update_, istr);
  feasst_deserialize(&is_block_, istr);
  feasst_deserialize_fstobj(&blocks_, istr);
}

void BiasTransitionMatrix::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(667, ostr);
  feasst_serialize_fstobj(ln_macro_prob_, ostr);
  feasst_serialize_fstobj(collection_, ostr);
  feasst_serialize(visits_, ostr);
  feasst_serialize(min_visits_, ostr);
  feasst_serialize(num_sweeps_, ostr);
  feasst_serialize(min_sweeps_, ostr);
  feasst_serialize(num_steps_to_update_, ostr);
  feasst_serialize(num_steps_since_update_, ostr);
  feasst_serialize(is_block_, ostr);
  feasst_serialize_fstobj(blocks_, ostr);
}

}  // namespace feasst
