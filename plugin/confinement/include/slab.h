
#ifndef FEASST_CONFINEMENT_SLAB_H_
#define FEASST_CONFINEMENT_SLAB_H_

#include "confinement/include/shape.h"
#include "core/include/arguments.h"

namespace feasst {

/**
  A slab is the intersection of two half spaces.
 */
class Slab : public Shape {
 public:
  Slab(
    /**
      dimension : Set the dimension which is finite.

      bound0 : Set a lower or upper value of the slab.

      bound1 : Set a lower or upper value of the slab.
     */
    const argtype &args = argtype());

  double nearest_distance(const Position& point) const override {
    return slab_->nearest_distance(point); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(485, ostr);
    feasst_serialize(slab_, ostr);
  }

  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Slab>(istr); }

  Slab(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(485 == version, version);
    // feasst_deserialize(slab_, istr);
    // HWH for unknown reasons, the above doesn't work
    int existing;
    istr >> existing;
    if (existing != 0) {
      slab_ = slab_->deserialize(istr);
    }
  }

  virtual ~Slab() {}

 private:
  const std::string class_name_ = "Slab";
  std::shared_ptr<Shape> slab_;
  Arguments args_;
};

inline std::shared_ptr<Slab> MakeSlab(
    const argtype& args = argtype()) {
  return std::make_shared<Slab>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_SLAB_H_
