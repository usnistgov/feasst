
#ifndef FEASST_CORE_UTILS_MATH_H_
#define FEASST_CORE_UTILS_MATH_H_

#include <vector>
#include <numeric>
#include <algorithm>

namespace feasst {

/** Return rounded double to nearest integer. This rounding is implemented
 *  as floor(x+0.5), such that feasstRound(-0.5) == 0. The cplusplus library
 *  round(-0.5) from math.h results in round(-0.5) == -1, such that rounding
 *  at the halfway point is away from zero. This breaks the first bin on
 *  histograms. */
int round(double x);

/// Return the average of a vector.
template<class T>
double average(const std::vector<T> &x) {
  return std::accumulate(x.begin(), x.end(), 0.)/static_cast<double>(x.size());
}

/// Return product of all elements of a vector.
template<class T>
T product(const std::vector<T> &vec) {
  if (static_cast<int>(vec.size()) == 0) {
    return 0;
  }
  T prod = 1;
  for (int i = 0; i < int(vec.size()); ++i) {
    prod *= vec[i];
  }
  return prod;
}

/// Return the minimum element of a vector.
template<class T>
T minimum(const std::vector<T> &vec) {
  return *std::min_element(vec.begin(), vec.begin() + vec.size());
}

}  // namespace feasst

#endif  // FEASST_CORE_UTILS_MATH_H_
