
#ifndef FEASST_CORE_UTILS_MATH_H_
#define FEASST_CORE_UTILS_MATH_H_

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include "core/include/constants.h"

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
    return static_cast<T>(0);
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
  return *std::min_element(vec.begin(), vec.end());
}

/// Return the minimum element of a 2D vector.
template<class T>
T minimum(const std::vector<std::vector<T> > &vec) {
  std::vector<T> mins;
  for (const std::vector<T>& vec1 : vec) {
    mins.push_back(*std::min_element(vec1.begin(), vec1.end()));
  }
  return *std::min_element(mins.begin(), mins.end());
}

/// Return the maximum element of a vector.
template<class T>
T maximum(const std::vector<T> &vec) {
  return *std::max_element(vec.begin(), vec.end());
}

/// Return the maximum element of a 2D vector.
template<class T>
T maximum(const std::vector<std::vector<T> > &vec) {
  std::vector<T> maxs;
  for (const std::vector<T>& vec1 : vec) {
    maxs.push_back(*std::max_element(vec1.begin(), vec1.end()));
  }
  return *std::max_element(maxs.begin(), maxs.end());
}

/// Compute the union of two vectors.
template<typename T>
std::vector<T> fst_union(const std::vector<T>& vec1,
    const std::vector<T>& vec2) {
  std::vector<T> both(vec1.size() + vec2.size());
  typename std::vector<T>::iterator iter;
  iter = std::set_union(vec1.begin(), vec1.end(),
                        vec2.begin(), vec2.end(),
                        both.begin());
  both.resize(iter - both.begin());
  return both;
}

/// Compute the difference of two vectors.
template<typename T>
std::vector<T> fst_difference(const std::vector<T>& vec1,
    const std::vector<T>& vec2) {
  std::vector<T> both(vec1.size() + vec2.size());
  typename std::vector<T>::iterator iter;
  iter = std::set_difference(vec1.begin(), vec1.end(),
                             vec2.begin(), vec2.end(),
                             both.begin());
  both.resize(iter - both.begin());
  return both;
}

/// Return the sign of the value.
/// Thanks to https://stackoverflow.com/questions/1903954/
/// is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/// Convert radians to degrees.
inline double radians_to_degrees(const double radians) {
  return radians/PI*180.;
}

/// Convert degrees to radians.
inline double degrees_to_radians(const double degrees) {
  return degrees/180.*PI;
}

/// Swap the values.
template <typename T>
inline void swap(T * val1, T * val2) {
  const double temp = *val1;
  *val1 = *val2;
  *val2 = temp;
}

/// Sort the values in order of increasing size.
template <typename T>
inline void sort(T * min, T * max) {
  if (*min > *max) {
    swap(min, max);
  }
}

/// Return the cumulative probability from a weighted series.
std::vector<double> cumulative_probability(const std::vector<double>& weights);

}  // namespace feasst

#endif  // FEASST_CORE_UTILS_MATH_H_
