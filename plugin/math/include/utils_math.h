
#ifndef FEASST_MATH_UTILS_MATH_H_
#define FEASST_MATH_UTILS_MATH_H_

#include <vector>
#include <numeric>  // accumulate
#include <algorithm>  // set_difference, set_union, max_element

namespace feasst {

// HWH rename utils_math -> utils

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
  for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
    prod *= vec[i];
  }
  return prod;
}

/// Return the minimum element of a vector.
template<class T>
T minimum(const std::vector<T> &vec) {
  // ASSERT(vec.size() > 0, "vector has no elements");
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

/// Return the sum of all elements in a vector.
template<class T>
T sum(const std::vector<T>& vec) {
  T sm = 0;
  for (const T& element : vec) {
    sm += element;
  }
  return sm;
}

/// Return the sum of all elements in a 2D vector.
template<class T>
T sum(const std::vector<std::vector<T> >& vec2) {
  T sm = 0;
  for (const std::vector<T>& vec : vec2) {
    sm += sum(vec);
  }
  return sm;
}

/// Return the sum of all elements in a 3D vector.
template<class T>
T sum(const std::vector<std::vector<std::vector<T> > >& vec3) {
  T sm = 0;
  for (const std::vector<std::vector<T> >& vec2 : vec3) {
    sm += sum(vec2);
  }
  return sm;
}

/// Return the sum of all elements in a 4D vector.
template<class T>
T sum(const std::vector<std::vector<std::vector<std::vector<T> > > >& vec4) {
  T sm = 0;
  for (const std::vector<std::vector<std::vector<T> > >& vec3 : vec4) {
    sm += sum(vec3);
  }
  return sm;
}

/** Return the minimum element of a 2D vector with smoothing.
    The returned indices must be global minimum between +/- num_smooth. */
template<class T>
std::vector<int> local_minimum_indices(const std::vector<T> &vec,
                                       const int num_smooth) {
  std::vector<int> mins;
  if (vec.size() == 0) {
    mins.push_back(0);
    return mins;
  }
  const int size = static_cast<int>(vec.size());
  for (int index = 1; index < size - 1; ++index) {
    int lower = index - num_smooth;
    if (lower < 0) lower = 0;
    int upper = index + num_smooth;
    if (upper > size - 1) upper = size - 1;
    auto min_iter = std::min_element(vec.begin() + lower,
                                     vec.begin() + upper + 1);
    const int minimum = std::distance(vec.begin(), min_iter);
    if (minimum == index) {
      mins.push_back(index);
    }
  }
  return mins;
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
std::vector<T> feasst_union(const std::vector<T>& vec1,
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
std::vector<T> feasst_difference(const std::vector<T>& vec1,
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
double radians_to_degrees(const double radians);

/// Convert degrees to radians.
double degrees_to_radians(const double degrees);

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

/// Return true if value is between bounds 1 and 2, inclusive.
bool is_in_interval(const double value,
                    double bound1,
                    double bound2);

/// Return the cumulative probability from a weighted series.
std::vector<double> cumulative_probability(const std::vector<double>& weights);

/// Return the volume in a spherical shell.
double spherical_shell_volume(const double lower,
  const double upper,
  const int dimension);

}  // namespace feasst

#endif  // FEASST_MATH_UTILS_MATH_H_
