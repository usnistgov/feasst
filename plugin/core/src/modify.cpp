
#include "core/include/debug.h"
#include "core/include/modify.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Modify> >& Modify::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Modify> >* ans =
     new std::map<std::string, std::shared_ptr<Modify> >();
  return *ans;
}

void Modify::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Modify> Modify::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Modify> Modify::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

class MapModifyUpdateOnly {
 public:
  MapModifyUpdateOnly() {
    ModifyUpdateOnly().deserialize_map()["ModifyUpdateOnly"] =
      std::make_shared<ModifyUpdateOnly>();
  }
};

static MapModifyUpdateOnly mapper_modify_update_only = MapModifyUpdateOnly();

std::shared_ptr<Modify> ModifyUpdateOnly::create(std::istream& istr) const {
  feasst_deserialize_version(istr);
  auto model = std::make_shared<ModifyUpdateOnly>();
  return model;
}

void ModifyUpdateOnly::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(1, ostr);
}

}  // namespace feasst
