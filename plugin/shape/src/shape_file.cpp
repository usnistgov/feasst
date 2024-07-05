#include <memory>
#include "utils/include/io.h"
#include "utils/include/file.h"
#include "utils/include/serialize.h"
#include "shape/include/sphere.h"
#include "shape/include/shape_file.h"
#include "shape/include/shape_union.h"
#include "shape/include/shape_intersect.h"

namespace feasst {

class MapShapeFile {
 public:
  MapShapeFile() {
    auto obj = std::make_shared<ShapeFile>();
    //auto obj = MakeShapeFile({{"file_name", install_dir() + "/plugin/shape/test/data/shape.txt"}});
    obj->deserialize_map()["ShapeFile"] = obj;
  }
};

static MapShapeFile mapper_shape_union_ = MapShapeFile();

ShapeFile::ShapeFile(argtype * args) : Shape() {
  class_name_ = "ShapeFile";
  DEBUG("args " << str(*args));
  std::string shape_file = str("shape_file", args, "");
  if (shape_file.empty()) {
    if (used("file_name", *args)) {
      WARN("ShapeFile::file_name was renamed to shape_file.");
      shape_file = str("file_name", args);
    } else {
      FATAL("ShapeFile::shape_file is a required argument.");
    }
  }
  std::ifstream file(shape_file);
  ASSERT(file.good(), "cannot find " << shape_file);
  std::string line;
  std::getline(file, line);
  DEBUG(line);
  std::pair<std::string, argtype> parsed_line = parse_line(line, NULL, NULL);
  DEBUG(parsed_line.first);
  DEBUG(str(parsed_line.second));
  if (line.empty()) return;

  // print every shape available in factory
  for (const auto& mp : deserialize_map()) {
    DEBUG(mp.first);
  }

  shape_ = this->factory(parsed_line.first, &parsed_line.second);
  while (std::getline(file, line)) {
    DEBUG(line);
    if (line.substr(0, 5) == "union") {
      DEBUG("union");
      line.replace(0, 6, "");
      DEBUG("new line: " << line);
      parsed_line = parse_line(line, NULL, NULL);
      shape_ = MakeShapeUnion(shape_,
        this->factory(parsed_line.first, &parsed_line.second));
    } else if (line.substr(0, 9) == "intersect") {
      DEBUG("intersect");
      line.replace(0, 10, "");
      DEBUG("new line: " << line);
      parsed_line = parse_line(line, NULL, NULL);
      shape_ = MakeShapeIntersect(shape_,
        this->factory(parsed_line.first, &parsed_line.second));
    }
  }
}
ShapeFile::ShapeFile(argtype args) : ShapeFile(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void ShapeFile::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_(ostr);
  feasst_serialize_version(6473, ostr);
  feasst_serialize_fstdr(shape_, ostr);
}

ShapeFile::ShapeFile(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(6473 == version, version);
  // HWH for unknown reasons, below template isn't working in this case.
  // feasst_deserialize_fstdr(shape, istr);
  int existing;
  istr >> existing;
  if (existing != 0) {
    shape_ = shape_->deserialize(istr);
  }
}

double ShapeFile::nearest_distance(const Position& point) const {
  return shape_->nearest_distance(point);
}

}  // namespace feasst
