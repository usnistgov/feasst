
#ifndef FEASST_UTILS_UTILS_H_
#define FEASST_UTILS_UTILS_H_

#include <cmath>
#include <vector>

namespace feasst {

/// Return true if value is found in list, and the index of that value in list.
template<class T>
bool find_in_list(const T value, const std::vector<T> &list,
  int * index  //!< last index in list where value was found
  ) {
  const int size = static_cast<int>(list.size());
  if (size == 0) return false;
  bool in = false;
  *index = -1;
  for (int i = 0; i < size; ++i) {
    if (list[i] == value) {
      in = true;
      *index = i;
    }
  }
  return in;
}

/// Same as above, but do not return the index.
template<class T>
bool find_in_list(const T value, const std::vector<T> &list) {
  int index;
  return find_in_list(value, list, &index);
}

/// Same as above, but use the first element of pair in list.
template<class T1, class T2>
bool find_in_list(const T1 value, const std::vector<std::pair<T1,T2> >& list,
  int * index  //!< last index in list where value was found
  ) {
  const int size = static_cast<int>(list.size());
  if (size == 0) return false;
  bool in = false;
  *index = -1;
  for (int i = 0; i < size; ++i) {
    if (list[i].first == value) {
      in = true;
      *index = i;
    }
  }
  return in;
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

/// Resize three dimensional vector.
/// Note to HWH: rewrite using variadic templates for arbitrary dimensions
template<class T>
void resize(const int xs, const int ys, const int zs,
  std::vector<std::vector<std::vector<T> > > *vec) {
  vec->resize(xs);
  for (int i = 0; i < static_cast<int>(vec->size()); ++i) {
    resize(ys, zs, &(*vec)[i]);
  }
}

/// Resize four dimensional vector.
template<class T>
void resize(const int dim0, const int dim1, const int dim2, const int dim3,
  std::vector<std::vector<std::vector<std::vector<T> > > > *vec) {
  vec->resize(dim0);
  for (int i = 0; i < static_cast<int>(vec->size()); ++i) {
    resize(dim1, dim2, dim3, &(*vec)[i]);
  }
}

/// Resize five dimensional vector.
template<class T>
void resize(const int dim0, const int dim1, const int dim2, const int dim3, const int dim4,
  std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > > *vec) {
  vec->resize(dim0);
  for (int i = 0; i < static_cast<int>(vec->size()); ++i) {
    resize(dim1, dim2, dim3, dim4, &(*vec)[i]);
  }
}

/// Fill two dimensional vector with a constant.
/// Note to HWH: rewrite using variadic templates for arbitrary dimensions
template<class T>
void fill(const T value, std::vector<std::vector<T> > *vec) {
  for (std::vector<T>& inner : *vec) {
    std::fill(inner.begin(), inner.end(), value);
  }
}

/// Fill three dimensional vector with a constant.
/// Note to HWH: rewrite using variadic templates for arbitrary dimensions
template<class T>
void fill(const T value, std::vector<std::vector<std::vector<T> > > *vec) {
  for (std::vector<std::vector<T> >& inner : *vec) {
    fill(value, &inner);
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

//template<class T>
//void feasst_reverse(std::vector<T> * vec) {
//  std::reverse(vec->begin(), vec->end());
//}
//
//template<class T>
//void feasst_reverse(std::vector<std::vector<T> > * vec) {
//  std::reverse(vec->begin(), vec->end());
//  for (std::vector<T>& vec1 : *vec) {
//    feasst_reverse(&vec1);
//  }
//}

/// Return true if the sorted vector contains a duplicate value.
template<class T>
bool has_duplicate(const std::vector<T>& vec) {
  std::vector<T> vec2 = vec;
  // HWH assume sorted, std::sort(vec2.begin(), vec2.end());
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
//    TRACE("size of vec1:" << vec1.size() << " != size of vec2:" << vec2.size());
    return false;
  }
  for (int index = 0; index < static_cast<int>(vec1.size()); ++index) {
    if (std::abs(vec1[index] - vec2[index]) > tolerance) {
//      TRACE(MAX_PRECISION << "vec1[" << index << "]:" << vec1[index] << " != "
//        << "vec2[" << index << "]:" << vec2[index]);
      return false;
    }
  }
  return true;
}

template<class T>
bool is_equal(const std::vector<std::vector<T> >& vec1,
              const std::vector<std::vector<T> >& vec2,
              const double tolerance = 1e-15) {
  if (vec1.size() != vec2.size()) return false;
  for (int index = 0; index < static_cast<int>(vec1.size()); ++index) {
    if (!is_equal(vec1[index], vec2[index], tolerance)) return false;
  }
  return true;
}

template<class T>
bool is_equal(const std::vector<std::vector<std::vector<T> > >& vec1,
              const std::vector<std::vector<std::vector<T> > >& vec2,
              const double tolerance = 1e-15) {
  if (vec1.size() != vec2.size()) return false;
  for (int index = 0; index < static_cast<int>(vec1.size()); ++index) {
    if (!is_equal(vec1[index], vec2[index])) return false;
  }
  return true;
}

template<class T>
bool is_equal(const std::vector<std::vector<std::vector<
              std::vector<T> > > >& vec1,
              const std::vector<std::vector<std::vector<
              std::vector<T> > > >& vec2,
              const double tolerance = 1e-15) {
  if (vec1.size() != vec2.size()) return false;
  for (int index = 0; index < static_cast<int>(vec1.size()); ++index) {
    if (!is_equal(vec1[index], vec2[index])) return false;
  }
  return true;
}

template<class T>
bool is_equal(const std::vector<std::vector<std::vector<std::vector<
              std::vector<T> > > > >& vec1,
              const std::vector<std::vector<std::vector<std::vector<
              std::vector<T> > > > >& vec2,
              const double tolerance = 1e-15) {
  if (vec1.size() != vec2.size()) return false;
  for (int index = 0; index < static_cast<int>(vec1.size()); ++index) {
    if (!is_equal(vec1[index], vec2[index])) return false;
  }
  return true;
}

}  // namespace feasst

#endif  // FEASST_UTILS_UTILS_H_
