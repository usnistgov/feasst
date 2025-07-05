
#ifndef FEASST_UTILS_INCLUDE_TIMER_RDTSC_H_
#define FEASST_UTILS_INCLUDE_TIMER_RDTSC_H_

#include <string>
#include <vector>
#include <cstdint>

namespace feasst {

/*
  TimerRDTSC is used to profile time spend on various parts of the code.
 */
class TimerRDTSC {
 public:
  /// Construct a timer with labels
  explicit TimerRDTSC(const int num = 0);

  /// Begin timing an index. If a previous index was timing, record previous.
  void start(
    /// If the index is -1, then simply end and record.
    const int index);

  /// Return the percentage of time spent in each index.
  std::vector<double> percents() const;

  /// Add an index.
  void add();

  /// Erase an index.
  void erase(const int index);

  void serialize(std::ostream& ostr) const;
  explicit TimerRDTSC(std::istream& istr);
  ~TimerRDTSC();

 private:
  std::vector<uint64_t> ticks_;
  int previous_index_ = -1;
  uint64_t previous_clock_ = 0;
};

}  // namespace feasst

#endif  // FEASST_UTILS_INCLUDE_TIMER_RDTSC_H_
