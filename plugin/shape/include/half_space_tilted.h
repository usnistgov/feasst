
#ifndef FEASST_SHAPE_HALF_SPACE_TILTED_H_
#define FEASST_SHAPE_HALF_SPACE_TILTED_H_

#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  A half space divides space by a plane (or line in 2D).
  This is a generalized implementation of HalfSpace.
  This plane is specified either by two points, or a normal vector and distance
  from origin.
  See https://mathworld.wolfram.com/Point-PlaneDistance.html
 */
class HalfSpaceTilted : public Shape {
 public:
  HalfSpaceTilted(
    /// A vector that is normal to the surface and points inside the half space.
    const Position& normal,
    /// The signed nearest distance of the surface to the origin.
    /// Negative if inside.
    const double distance_from_origin);

  /**
    Alternatively, construct the planar surface with two points.
    The first is on the plane, and the second is inside the half space.
    The vector connecting these two points is perpendicular to the plane.
   */
  HalfSpaceTilted(const Position& point0, const Position& point1);

  /// Return the unit normal vector.
  const Position& unit_normal() const { return unit_normal_; }

  /// Return the distance from the origin.
  const double distance_from_origin() const { return distance_from_origin_; }

  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<HalfSpaceTilted>(istr); }
  explicit HalfSpaceTilted(std::istream& istr);
  virtual ~HalfSpaceTilted() {}

 protected:
  void serialize_half_space_tilted_(std::ostream& ostr) const;

 private:
  Position unit_normal_;
  double distance_from_origin_;
  void init_(const Position& unit_normal, const double distance_from_origin);
};

inline std::shared_ptr<HalfSpaceTilted> MakeHalfSpaceTilted(
    const Position& point0,
    const Position& point1) {
  return std::make_shared<HalfSpaceTilted>(point0, point1); }

inline std::shared_ptr<HalfSpaceTilted> MakeHalfSpaceTilted(
    const Position& normal,
    const double distance_from_origin) {
  return std::make_shared<HalfSpaceTilted>(normal, distance_from_origin); }

}  // namespace feasst

#endif  // FEASST_SHAPE_HALF_SPACE_TILTED_H_
