#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "steppers/include/seek_modify.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

std::vector<int> SeekModify::index(const std::string class_name,
    const MonteCarlo& mc) const {
  DEBUG("Looking for " << class_name);
  DEBUG("num " << mc.num_modifiers());
  for (int index1 = 0; index1 < mc.num_modifiers(); ++index1) {
    const Modify& mod = mc.modify(index1);
    std::string name = mod.class_name();
    DEBUG("found " << name);
    if (name == class_name) {
      return {index1, -1};
    }
    if (name == "ModifyFactory") {
      const auto& mods = mod.modifiers();
      for (int index2 = 0; index2 < static_cast<int>(mods.size()); ++index2) {
        name = mods[index2]->class_name();
        if (name == class_name) {
          return {index1, index2};
        }
      }
    }
  }
  return {-1, -1};
}

}  // namespace feasst
