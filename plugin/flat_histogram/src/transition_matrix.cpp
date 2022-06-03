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
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/transition_matrix.h"

namespace feasst {

TransitionMatrix::TransitionMatrix(argtype args) : TransitionMatrix(&args) {
  FEASST_CHECK_ALL_USED(args);
}
TransitionMatrix::TransitionMatrix(argtype * args) {
  class_name_ = "TransitionMatrix";
  min_visits_ = integer("min_visits", args, 100);
  min_sweeps_ = integer("min_sweeps", args);
  reset_sweeps_ = integer("reset_sweeps", args, -1);
  collection_ = CollectionMatrix(args);
  average_visits_ = integer("average_visits", args, 0);
  new_sweep_ = integer("new_sweep", args, 0);
}

void TransitionMatrix::update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const Macrostate& macro) {
  DEBUG("macro old/new " << macrostate_old << " " << macrostate_new);
  DEBUG("is_accepted " << is_accepted);
  DEBUG("is_endpoint " << is_endpoint);
  const int index = macrostate_new - macrostate_old + 1;
  ASSERT(index >= 0 and index <= 2, "index(" << index << ") must be 0, 1 or 2");
  if (is_accepted && (macrostate_old != macrostate_new)) {
    int vindex = 0;
    if (macrostate_old < macrostate_new) {
      vindex = 1;
    }
    ++visits_[macrostate_new][vindex];
    // double count endpoints because transitions beyond bounds are impossible
    if (is_endpoint) {
      ++visits_[macrostate_new][1-vindex];
    }
    if (new_sweep_ != 0) {
      num_sweeps_ = min_vis_calc_(macro);
    }
  }
  double metropolis_prob = std::min(1., std::exp(ln_metropolis_prob));
  DEBUG("macrostate_old " << macrostate_old << " index " << index);
  DEBUG("metropolis_prob " << metropolis_prob);
  if (index == 0) {
    collection_.increment(macrostate_old, 0, metropolis_prob);
    collection_.increment(macrostate_old, 1, 0.);
  } else if (index == 1) {
    collection_.increment(macrostate_old, 0, 0.);
    collection_.increment(macrostate_old, 1, 0.);
  } else if (index == 2) {
    collection_.increment(macrostate_old, 0, 0.);
    collection_.increment(macrostate_old, 1, metropolis_prob);
  } else {
    FATAL("unrecognized index: " << index);
  }
}

void TransitionMatrix::resize(const int size) {
  ln_prob_.resize(size);
  feasst::resize(size, 2, &visits_);
  collection_.resize(size);
}

std::string TransitionMatrix::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << "\"num_sweeps\":" << num_sweeps_ << ",";
  //TRACE("matrix," << feasst_str(collection_.matrix()));
  return ss.str();
}

std::string TransitionMatrix::write_per_bin_header() const {
  std::stringstream ss;
  ss << Bias::write_per_bin_header() << ",";
  ss << "visits0,visits1,";
  ss << collection_.write_per_bin_header();
  return ss.str();
}

std::string TransitionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << ",";
  ss << MAX_PRECISION << visits_[bin][0] << "," << visits_[bin][1] << ",";
  ss << collection_.write_per_bin(bin);
  return ss.str();
}

int TransitionMatrix::min_vis_calc_(const Macrostate& macro) const {
  int min_vis = 1e9;
  for (int bin = macro.soft_min(); bin < macro.soft_max() + 1; ++bin) {
    int vis = *std::min_element(visits_[bin].begin(), visits_[bin].end());
    if (vis < min_vis) {
      min_vis = vis;
    }
  }
  // ASSERT(min_vis != 1e9, "error");
  return min_vis;
}

void TransitionMatrix::infrequent_update(const Macrostate& macro) {
  DEBUG("TransitionMatrix::infrequent_update_()");
  DEBUG("update the macrostate distribution");
  collection_.compute_ln_prob(&ln_prob_);

  DEBUG("update the number of sweeps");
  if (new_sweep_ == 0) {
    if (min_vis_calc_(macro) >= min_visits_) {
      double sum_vis = 0;
      for (int bin = macro.soft_min(); bin < macro.soft_max() + 1; ++bin) {
        sum_vis += visits_[bin][0] + visits_[bin][1];
      }
      if (sum_vis >= average_visits_*(macro.soft_max() - macro.soft_min() + 1)) {
        ++num_sweeps_;
        feasst::fill(0, &visits_);
        //std::fill(visits_.begin(), visits_.end(), 0);
      }
    }
  }
  if (num_sweeps_ == reset_sweeps_) {
     increment_phase();
  }

  DEBUG("check if complete");
  if (num_sweeps_ >= min_sweeps_) {
    set_complete_();
  } else {
    set_incomplete_();
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
  ASSERT(667 == version || version == 668 || version == 669,
    "mismatch version: " << version);
  feasst_deserialize_fstobj(&ln_prob_, istr);
  feasst_deserialize_fstobj(&collection_, istr);
  feasst_deserialize(&visits_, istr);
  feasst_deserialize(&min_visits_, istr);
  feasst_deserialize(&num_sweeps_, istr);
  feasst_deserialize(&min_sweeps_, istr);
  feasst_deserialize(&reset_sweeps_, istr);
  if (version > 668) {
    feasst_deserialize(&average_visits_, istr);
    feasst_deserialize(&new_sweep_, istr);
  }
//  if (version == 667) {  // HWH backwards compatability
//    int num_blocks, iter_block;
//    bool is_block;
//    std::vector<TransitionMatrix> blocks;
//    feasst_deserialize(&num_blocks, istr);
//    feasst_deserialize(&is_block, istr);
//    feasst_deserialize(&iter_block, istr);
//    feasst_deserialize_fstobj(&blocks, istr);
//  }
}

void TransitionMatrix::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(669, ostr);
  feasst_serialize_fstobj(ln_prob_, ostr);
  feasst_serialize_fstobj(collection_, ostr);
  feasst_serialize(visits_, ostr);
  feasst_serialize(min_visits_, ostr);
  feasst_serialize(num_sweeps_, ostr);
  feasst_serialize(min_sweeps_, ostr);
  feasst_serialize(reset_sweeps_, ostr);
  feasst_serialize(average_visits_, ostr);
  feasst_serialize(new_sweep_, ostr);
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
    INFO("visits not equal: " << feasst_str(visits_) << " vs "
      << feasst_str(transition_matrix.visits_));
    return false;
  }
  if (min_visits_ != transition_matrix.min_visits_) return false;
  if (num_sweeps_ != transition_matrix.num_sweeps_) return false;
  if (min_sweeps_ != transition_matrix.min_sweeps_) return false;
  if (reset_sweeps_ != transition_matrix.reset_sweeps_) return false;
  if (average_visits_ != transition_matrix.average_visits_) return false;
  if (new_sweep_ != transition_matrix.new_sweep_) return false;
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

void TransitionMatrix::set_cm(const int macro, const Bias& bias) {
  collection_.set(macro, bias.cm().matrix()[macro]);
  visits_[macro][0] = bias.visits(macro, 0);
  visits_[macro][1] = bias.visits(macro, 1);
}

void TransitionMatrix::set_cm(const CollectionMatrix& cm) {
  collection_ = cm;
  const int size = collection_.matrix().size();
  ln_prob_.resize(size);
  feasst::resize(size, 2, &visits_);
  collection_.compute_ln_prob(&ln_prob_);
}

int TransitionMatrix::num_iterations(const int state) const {
  if (new_sweep_ == 0 || state == -1) {
    return num_sweeps_;
  }
  return std::min(visits_[state][0], visits_[state][1]);
}

}  // namespace feasst
