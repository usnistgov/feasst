#include <algorithm>
#include "utils/include/serialize.h"
#include "morph/include/macrostate_morph.h"

namespace feasst {

MacrostateMorph::MacrostateMorph(
    const std::vector<std::vector<int> > grow_seq,
    const Histogram& histogram,
    const argtype& args) : Macrostate(histogram, args) {
  grow_seq_ = grow_seq;
  ASSERT(grow_seq.size() > 1, "The particle growth sequence is of size: " <<
    grow_seq.size() << " but must be 2 or more");
  ASSERT(grow_seq[0].size() > 0, "size of grow_seq[0]: " << grow_seq[0].size());
  for (int type : grow_seq.front()) {
    num_first_.push_back(ConstrainNumParticles({{"type", feasst::str(type)}}));
  }
}

bool MacrostateMorph::is_row_all_populated_(const int row,
    const Matrix& matrix) const {
  for (const double val : matrix.matrix()[row]) {
    if (val < 1e-8) {
      return false;
    }
  }
  return true;
}

double MacrostateMorph::value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  // begin by initializing the number matrix
  // number of particles of type, where row follows macrostate (steps),
  // and column is for each morphing particle.
  // HWH optimize: num_matrix_ could be private data, but method cant be const
  Matrix num_mat_;
  num_mat_.set_size(grow_seq_.size(), grow_seq_[0].size());
  num_mat_.zero();
  for (int col = 0; col < num_mat_.num_columns(); ++col) {
    num_mat_.set_value(0, col,
      num_first_[col].num_particles(system, acceptance));
//    num_mat_.set_value(num_mat_.num_rows() - 1, col,
//      system.configuration().num_particles_of_type(grow_seq_.back()[col]));
  }
  for (int row = 1; row < num_mat_.num_rows(); ++row) {
    for (int col = 0; col < num_mat_.num_columns(); ++col) {
      num_mat_.set_value(row, col,
        system.configuration().num_particles_of_type(grow_seq_[row][col]));
    }
  }
  int state = -1;
  if (is_row_all_populated_(0, num_mat_)) {
    state = 1;
  } else {
//      ASSERT(num_first == 0, "There are num_first: "
//        << num_first << " partial particles");
    for (int row = 1; row < num_mat_.num_rows() - 1; ++row) {
      if (is_row_all_populated_(row, num_mat_)) {
        state = row + 1;
        break;
      }
    }
    if (state == -1) {
      state = 0;
    }
  }
  DEBUG("state: " << state);
  DEBUG("num_mat: " << num_mat_.str());
  const double val = num_mat_.value(num_mat_.num_rows()-1, 0) +
    static_cast<double>(state)/
    static_cast<double>(grow_seq_.size());
  DEBUG("val: " << val);
  return val;
}

class MapMacrostateMorph {
 public:
  MapMacrostateMorph() {
    auto obj = MakeMacrostateMorph({{1}, {0}},
      Histogram({{"width", "1"}, {"max", "1"}}));
    obj->deserialize_map()["MacrostateMorph"] = obj;
  }
};

static MapMacrostateMorph mapper_ = MapMacrostateMorph();

std::shared_ptr<Macrostate> MacrostateMorph::create(std::istream& istr) const {
  return std::make_shared<MacrostateMorph>(istr);
}

MacrostateMorph::MacrostateMorph(std::istream& istr)
  : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5337, "version mismatch: " << version);
  feasst_deserialize(&grow_seq_, istr);
  feasst_deserialize_fstobj(&num_first_, istr);
}

void MacrostateMorph::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(5337, ostr);
  feasst_serialize(grow_seq_, ostr);
  feasst_serialize_fstobj(num_first_, ostr);
}

}  // namespace feasst
