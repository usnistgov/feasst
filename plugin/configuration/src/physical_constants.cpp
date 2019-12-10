
#include "utils/include/serialize.h"
#include "configuration/include/physical_constants.h"

namespace feasst {

std::map<std::string, std::shared_ptr<PhysicalConstants> >& PhysicalConstants::deserialize_map() {
  static std::map<std::string, std::shared_ptr<PhysicalConstants> >* ans =
     new std::map<std::string, std::shared_ptr<PhysicalConstants> >();
  return *ans;
}

std::shared_ptr<PhysicalConstants> PhysicalConstants::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

class MapCODATA2018 {
 public:
  MapCODATA2018() {
    auto obj = MakeCODATA2018();
    obj->deserialize_map()["CODATA2018"] = obj;
  }
};

static MapCODATA2018 mapper_codata2018_ = MapCODATA2018();

class MapCODATA2014 {
 public:
  MapCODATA2014() {
    auto obj = MakeCODATA2014();
    obj->deserialize_map()["CODATA2014"] = obj;
  }
};

static MapCODATA2014 mapper_codata2014_ = MapCODATA2014();

class MapCODATA2010 {
 public:
  MapCODATA2010() {
    auto obj = MakeCODATA2010();
    obj->deserialize_map()["CODATA2010"] = obj;
  }
};

static MapCODATA2010 mapper_codata2010_ = MapCODATA2010();

} // namespace feasst
