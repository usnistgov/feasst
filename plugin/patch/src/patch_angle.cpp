#include <cmath>
#include "utils/include/serialize.h"
#include "patch/include/patch_angle.h"
#include "math/include/utils_math.h"

namespace feasst {

class MapPatchAngle {
 public:
  MapPatchAngle() {
    auto obj = std::make_shared<PatchAngle>();
    obj->deserialize_map()["patch_angle"] = obj;
  }
};

static MapPatchAngle mapper_patch_angle_ = MapPatchAngle();

void PatchAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2957, ostr);
}

PatchAngle::PatchAngle(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2957, "mismatch version: " << version);
}

class MapDirector {
 public:
  MapDirector() {
    auto obj = std::make_shared<Director>();
    obj->deserialize_map()["director"] = obj;
  }
};

static MapDirector mapper_director_ = MapDirector();

void Director::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(1948, ostr);
}

Director::Director(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1948, "mismatch version: " << version);
}

double Director::compute(const int type, const ModelParams& model_params) {
  return model_params.select("director").value(type);
}

class MapCosPatchAngle {
 public:
  MapCosPatchAngle() {
    auto obj = std::make_shared<CosPatchAngle>();
    obj->deserialize_map()["cos_patch_angle"] = obj;
  }
};

static MapCosPatchAngle mapper_cos_patch_angle_ = MapCosPatchAngle();

void CosPatchAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(3967, ostr);
}

CosPatchAngle::CosPatchAngle(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3967, "mismatch version: " << version);
}

double CosPatchAngle::compute(const int type, const ModelParams& model_params) {
  const double angle = model_params.select("patch_angle").value(type);
  return std::cos(degrees_to_radians(angle));
}

class MapSpherocylinderLength {
 public:
  MapSpherocylinderLength() {
    auto obj = std::make_shared<SpherocylinderLength>();
    obj->deserialize_map()["spherocylinder_length"] = obj;
  }
};

static MapSpherocylinderLength mapper_length_ = MapSpherocylinderLength();

void SpherocylinderLength::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2073, ostr);
}

SpherocylinderLength::SpherocylinderLength(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2073, "mismatch version: " << version);
}

double SpherocylinderLength::compute(const int type, const ModelParams& model_params) {
  return model_params.select("spherocylinder_length").value(type);
}


}  // namespace feasst
