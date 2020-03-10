
#ifndef FEASST_UTILS_INCLUDE_TIMER_H_
#define FEASST_UTILS_INCLUDE_TIMER_H_

#include <time.h>
#include <string>
#include <vector>
#include <map>

namespace feasst {

inline double cpu_hours() {
  return static_cast<double>(clock())
         /static_cast<double>(CLOCKS_PER_SEC)/60./60.;
}

inline double cpu_hours(const clock_t& clock) {
  return static_cast<double>(clock)
         /static_cast<double>(CLOCKS_PER_SEC)/60./60.;
}

/*
  Timer is used to profile time spend on various parts of the code.
 */
class Timer {
 public:
  /// Add a new profile name. Return the index.
  int add(const std::string& name);

  /// Set the profile to the index and record the time.
  /// When a new index is started, or timer is ended, then the elapsed time
  /// is recorded.
  void start(const int index);

  /// End the current timer. Record previous index if exists.
  void end();

  /// Return the sum of clock ticks.
  clock_t sum_clocks() const;

  /// Return the number of clock ticks for an index.
  clock_t clocks(const int index) const { return clocks_[index]; }

  /// Return the total CPU hours of the index.
  double hours(const int index) const;

  /// Return the total CPU hours of the name.
  double hours(const std::string name) const {
    return hours(name_to_index_(name)); }

  /// Return the total CPU hours of all indices.
  double hours() const;

  std::string str() const;

  /// The total clocks in this timer for name should be the same as total clocks
  /// in timer for all names. Return the missing percent.
  double missing_percent(const std::string name, const Timer& timer) const;

  /// Serialize object.
//  void serialize(std::ostream& ostr) const;

  /// Deserialize object.
//  explicit Timer(std::istream& istr);

 private:
  std::vector<clock_t> clocks_;
  std::vector<std::string> names_;
//  std::map<std::string, clock_t> profiles_;
//  std::string previous_name_ = "";
  int previous_index_ = -1;
  clock_t previous_clock_ = 0;

  int name_to_index_(const std::string name) const;
};

}  // namespace feasst

#endif  // FEASST_UTILS_INCLUDE_TIMER_H_
