
#ifndef FEASST_UTILS_UTILS_H_
#define FEASST_UTILS_UTILS_H_

#include <vector>
#include <numeric>
#include <algorithm>

namespace feasst {

/// Return if value is found in list, and the index of that value in list.
template<class T>
bool find_in_list(const T value, const std::vector<T> &list,
  int * index  //!< last index in list where value was found
  ) {
  if (list.size() == 0) {
    return false;
  }
  bool in = false;
  *index = -1;
  for (int i = 0; i < static_cast<int>(list.size()); ++i) {
    if (list[i] == value) {
      in = true;
      *index = i;
    }
  }
  return in;
}

/// Return if value is found in list
template<class T>
bool find_in_list(const T value, const std::vector<T> &list) {
  int index;
  return find_in_list(value, list, &index);
}

/// Resize two dimensional vector to size xs, ys.
/// Note to HWH: rewrite using variadic templates for arbitrary dimensions
template<class T>
void resize(const int xs, const int ys, std::vector<std::vector<T> > *vec) {
  vec->resize(xs);
  for (int i = 0; i < static_cast<int>(vec->size()); ++i) {
    (*vec)[i].resize(ys);
  }
}

/// Return the total number of elements in a multidimensional vector.
template<class T>
int num_elements(const std::vector<std::vector<T> > vec) {
  int num = 0;
  for (const std::vector<T>& element : vec) {
    num += static_cast<int>(element.size());
  }
  return num;
}

template<class T>
void feasst_reverse(std::vector<T> * vec) {
  std::reverse(vec->begin(), vec->end());
}

template<class T>
void feasst_reverse(std::vector<std::vector<T> > * vec) {
  std::reverse(vec->begin(), vec->end());
  for (std::vector<T>& vec1 : *vec) {
    feasst_reverse(&vec1);
  }
}

template<class T>
bool has_duplicate(const std::vector<T>& vec) {
  std::vector<T> vec2 = vec;
  std::sort(vec2.begin(), vec2.end());
  for (int i = 1; i < static_cast<int>(vec2.size()); ++i) {
    if (vec2[i] == vec2[i - 1]) {
      return true;
    }
  }
  return false;
}

template<class T>
bool is_equal(const std::vector<T>& vec1,
              const std::vector<T>& vec2,
              const double tolerance = 1e-15) {
  if (vec1.size() != vec2.size()) {
    // DEBUG("size of vec1:" << vec1.size() << " != size of vec2:" <<
    //       vec2.size());
    return false;
  }
  for (int i = 1; i < static_cast<int>(vec1.size()); ++i) {
    if (std::abs(vec1[i] - vec2[i]) > tolerance) {
      // DEBUG(MAX_PRECISION << "vec1[" << i << "]:" << vec1[i] << " != " <<
      //       "vec2[" << i << "]:" << vec2[i]);
      return false;
    }
  }
  return true;
}

}  // namespace feasst

#endif  // FEASST_UTILS_UTILS_H_
