
#ifndef FEASST_SYSTEM_SYNCHRONIZE_DATA_H_
#define FEASST_SYSTEM_SYNCHRONIZE_DATA_H_

#include <cstdint>
#include <utility>
#include <vector>
#include <ostream>

namespace feasst {

typedef std::vector<std::vector<std::vector<std::vector<double> > > > vec4;
typedef std::vector<vec4> vec5;
typedef std::vector<vec5> vec6;
typedef std::vector<std::pair<int, std::vector<double> > > vpv;
typedef std::vector<std::pair<int, vpv> > vpvpv;
typedef std::vector<std::pair<int, vpvpv> > vpvpvpv;
typedef std::vector<std::pair<int, vpvpvpv> > vpvpvpvpv;
typedef std::vector<vpvpv> vvpvpv;
typedef std::vector<vvpvpv> vvvpvpv;

/**
  This is a generic data container used to synchronized derived classes from
  base classes.
 */
class SynchronizeData {
 public:
  SynchronizeData() {}

  /// Return 1D double precision data.
  const std::vector<double>& dble_1D() const { return dble_1D_; }

  /// Return 1D 32-bit integer data.
  const std::vector<int>& int_1D() const { return int_1D_; }

  /// Return 1D 64-bit integer data.
  const std::vector<int64_t>& int64_1D() const { return int64_1D_; }

  /// Get 1D double precision data.
  std::vector<double> * get_dble_1D() { return &dble_1D_; }

  /// Get 1D 32-bit integer data.
  std::vector<int> * get_int_1D() { return &int_1D_; }

  /// Return 2D 32-bit integer data.
  const std::vector<std::vector<int> >& int_2D() const { return int_2D_; }

  /// Get 1D 32-bit integer data.
  std::vector<std::vector<int> > * get_int_2D() { return &int_2D_; }

  /// Get 1D 64-bit integer data.
  std::vector<int64_t> * get_int64_1D() { return &int64_1D_; }

  /// Return 2D data.
  const std::vector<std::vector<double> >& dble_2D() const { return dble_2D_; }

  /// Get 2D data.
  std::vector<std::vector<double> > * get_dble_2D() { return &dble_2D_; }

  /// Return 3D data.
  const std::vector<std::vector<std::vector<double> > >& dble_3D() const {
    return dble_3D_; }

  /// Get 3D data.
  std::vector<std::vector<std::vector<double> > > * get_dble_3D() {
    return &dble_3D_; }

  /// Return 4D data.
  const vec4& dble_4D() const { return dble_4D_; }

  /// Get 4D data.
  vec4 * get_dble_4D() { return &dble_4D_; }

  /// Return 5D data.
  const vec5& dble_5D() const { return dble_5D_; }

  /// Get 5D data.
  vec5 * get_dble_5D() { return &dble_5D_; }

  /// Return 6D data.
  const vec6& dble_6D() const { return dble_6D_; }

  /// Get 6D data.
  vec6 * get_dble_6D() { return &dble_6D_; }

  /// Return vpvpvpvpv data.
  const vpvpvpvpv& get_const_vpvpvpvpv() const { return vpvpvpvpv_; }

  /// Get vpvpvpvpv data.
  vpvpvpvpv * get_vpvpvpvpv() { return &vpvpvpvpv_; }

  /// Return vvvpvpv data.
  const vvvpvpv& get_const_vvvpvpv() const { return vvvpvpv_; }

  /// Get vvvpvpv data.
  vvvpvpv * get_vvvpvpv() { return &vvvpvpv_; }

  void serialize(std::ostream& ostr) const;
  explicit SynchronizeData(std::istream& istr);

 private:
  std::vector<double> dble_1D_;
  std::vector<int> int_1D_;
  std::vector<int64_t> int64_1D_;
  std::vector<std::vector<int> > int_2D_;
  std::vector<std::vector<double> > dble_2D_;
  std::vector<std::vector<std::vector<double> > > dble_3D_;
  vec4 dble_4D_;
  vec5 dble_5D_;
  vec6 dble_6D_;
  vpvpvpvpv vpvpvpvpv_;
  vvvpvpv vvvpvpv_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_SYNCHRONIZE_DATA_H_
