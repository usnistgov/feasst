#include <sstream>
#include <iostream>
#include <numeric>
#include "utils/include/timer.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"

namespace feasst {

int Timer::add(const std::string& name) {
  names_.push_back(name);
  clocks_.push_back(0);
  return static_cast<int>(clocks_.size() - 1);
}

void Timer::start(const int index) {
  const clock_t current_clock = clock();
  if (previous_index_ != -1) {
    ASSERT(previous_clock_ != 0, "uninitialized");
    clocks_[previous_index_] += current_clock - previous_clock_;
  }
  previous_index_ = index;
  previous_clock_ = current_clock;
}

void Timer::end() {
  if (previous_index_ != -1) {
    const clock_t current_clock = clock();
    ASSERT(previous_clock_ != 0, "uninitialized");
    clocks_[previous_index_] += current_clock - previous_clock_;
  }
  previous_index_ = -1;
  previous_clock_ = 0;
}

clock_t Timer::sum_clocks() const {
  // HWH debug: does "0" on std::accumulate work OK?
  return std::accumulate(clocks_.begin(), clocks_.end(), 0);
}

double Timer::hours(const int index) const {
  return cpu_hours(clocks_[index]);
}

double Timer::hours() const {
  return cpu_hours(sum_clocks());
}

std::string Timer::str() const {
  std::stringstream ss;
  ss << "name percent clocks hours" << std::endl;
  for (int i = 0; i < static_cast<int>(clocks_.size()); ++i) {
    ss << names_[i] << "\t"
       << 100*hours(i)/hours() << "\t"
       << clocks_[i] << "\t"
       << hours(i) << "\t"
       << std::endl;
  }
  return ss.str();
}

int Timer::name_to_index_(const std::string name) const {
  int index = -1;
  const bool found = find_in_list(name, names_, &index);
  ASSERT(found, "name " << name << " not found.");
  return index;
}

double Timer::missing_percent(const std::string name,
                              const Timer& timer) const {
  const double major = hours(name);
  const double percent = 100.*(major - timer.hours())/major;
  // ASSERT(percent >= 0, "percent: " << percent << " should be positive.");
  return percent;
}

}  // namespace feasst
