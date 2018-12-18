
#ifndef FEASST_CORE_UTILS_H_
#define FEASST_CORE_UTILS_H_

#include <vector>
#include <numeric>

namespace feasst {

/// Return if value is found in list, and the index of that value in list.
template<class T>
bool find_in_list(const T value, const std::vector<T> &list,
  int * index  //!< last index in list where value was found
  ) {
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

}  // namespace feasst

#endif  // FEASST_CORE_UTILS_H_
