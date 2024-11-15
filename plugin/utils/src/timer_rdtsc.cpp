#include <x86intrin.h>
#include <numeric>  // accumulate
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/timer_rdtsc.h"

namespace feasst {

TimerRDTSC::TimerRDTSC(const int num) {
  ticks_.resize(num);
}
TimerRDTSC::~TimerRDTSC() {}

void TimerRDTSC::start(const int index) {
  const uint64_t current_clock = __rdtsc();
  if (previous_index_ != -1) {
    ASSERT(previous_clock_ != 0, "uninitialized");
    ticks_[previous_index_] += current_clock - previous_clock_;
  }
  previous_index_ = index;
  previous_clock_ = current_clock;
}

std::vector<double> TimerRDTSC::percents() const {
  const uint64_t total = std::accumulate(ticks_.begin(), ticks_.end(), 0.);
  std::vector<double> pers;
  for (const uint64_t val : ticks_) {
    pers.push_back(static_cast<double>(val)/total*100.);  
  }
  return pers;
}

void TimerRDTSC::add() {
  ticks_.push_back(0);
}

void TimerRDTSC::erase(const int index) {
  if (index == previous_index_) {
    start(-1);
  }
  ticks_.erase(ticks_.begin() + index);
}

void TimerRDTSC::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2346, ostr);
  feasst_serialize(static_cast<int>(ticks_.size()), ostr);
}

TimerRDTSC::TimerRDTSC(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2346, "version mismatch: " << version);
  int size;
  feasst_deserialize(&size, istr);
  ticks_.resize(size);
}

}  // namespace feasst
