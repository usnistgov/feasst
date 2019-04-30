
#include "flat_histogram/include/macrostate_accumulator.h"

namespace feasst {

void MacrostateAccumulatorFactory::resize(const int size) {
  trackers_.clear();
  for (int index = 0; index < size; ++index) {
    std::vector<std::shared_ptr<MacrostateAccumulator> > trk;
    for (std::shared_ptr<MacrostateAccumulator> type : tracker_types_) {
      trk.push_back(type->deep_copy());
    }
    trackers_.push_back(trk);
  }
}

void MacrostateAccumulatorFactory::update(
    const int bin,
    const System * system,
    const Criteria * criteria) {
  for (std::shared_ptr<MacrostateAccumulator> tracker : trackers_[bin]) {
    tracker->update(system, criteria);
  }
}

std::string MacrostateAccumulatorFactory::write_per_bin_header() const {
  std::stringstream ss;
  for (auto tracker : trackers_[0]) {
    ss << tracker->name() << " "
       << tracker->name() << "_stdev ";
  }
  return ss.str();
}

std::string MacrostateAccumulatorFactory::write_per_bin(const int bin) const {
  std::stringstream ss;
  for (auto tracker : trackers_[bin]) {
    ss << tracker->accumulator().average() << " "
       << tracker->accumulator().block_stdev();
  }
  return ss.str();
}

}  // namespace feasst
