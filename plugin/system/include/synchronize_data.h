
#ifndef FEASST_SYSTEM_SYNCHRONIZE_DATA_H_
#define FEASST_SYSTEM_SYNCHRONIZE_DATA_H_

#include <vector>
#include <sstream>

namespace feasst {

typedef std::vector<std::vector<std::vector<std::vector<
  std::vector<double> > > > > vec5;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<
  std::vector<double> > > > > > vec6;

/**
  This is a generic data container used to synchronized derived classes from
  base classes.
 */
class SynchronizeData {
 public:
  SynchronizeData() {}

  /// Return 1D data.
  const std::vector<double>& dble_1D() const { return dble_1D_; }

  /// Return 1D data.
  const std::vector<int64_t>& int64_1D() const { return int64_1D_; }

  /// Return 2D data.
  const std::vector<std::vector<double> >& dble_2D() const { return dble_2D_; }

  /// Return 5D data.
  const vec5& dble_5D() const { return dble_5D_; }

  /// Return 6D data.
  const vec6& dble_6D() const { return dble_6D_; }

  /// Get 1D data.
  std::vector<double> * get_dble_1D() { return &dble_1D_; }

  /// Get 1D data.
  std::vector<int64_t> * get_int64_1D() { return &int64_1D_; }

  /// Get 2D data.
  std::vector<std::vector<double> > * get_dble_2D() { return &dble_2D_; }

  /// Get 5D data.
  vec5 * get_dble_5D() { return &dble_5D_; }

  /// Get 6D data.
  vec6 * get_dble_6D() { return &dble_6D_; }

  void serialize(std::ostream& ostr) const;
  explicit SynchronizeData(std::istream& istr);

 private:
  std::vector<double> dble_1D_;
  std::vector<int64_t> int64_1D_;
  std::vector<std::vector<double> > dble_2D_;
  vec5 dble_5D_;
  vec6 dble_6D_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_SYNCHRONIZE_DATA_H_
