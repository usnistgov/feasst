
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/convert_to_ref_potential.h"

namespace feasst {

ConvertToRefPotential::ConvertToRefPotential(argtype * args) {
  class_name_ = "ConvertToRefPotential";
  potential_index_ = integer("potential_index", args, 0);
  cutoff_ = dble("cutoff", args, -1);
  use_cell_ = boolean("use_cell", args, false);
}
ConvertToRefPotential::ConvertToRefPotential(argtype args) : ConvertToRefPotential(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapConvertToRefPotential {
 public:
  MapConvertToRefPotential() {
    auto obj = MakeConvertToRefPotential();
    obj->deserialize_map()["ConvertToRefPotential"] = obj;
  }
};

static MapConvertToRefPotential mapper_ConvertToRefPotential = MapConvertToRefPotential();

ConvertToRefPotential::ConvertToRefPotential(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8473, "mismatch version: " << version);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&cutoff_, istr);
  feasst_deserialize(&use_cell_, istr);
}

void ConvertToRefPotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(8473, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(cutoff_, ostr);
  feasst_serialize(use_cell_, ostr);
}

void ConvertToRefPotential::run(MonteCarlo * mc) {
  const Potential& pot = mc->system().potential(potential_index_);
  std::stringstream ss;
  pot.serialize(ss);
  std::shared_ptr<Potential> ref = std::make_shared<Potential>(ss);
  if (cutoff_ > 0) {
    const Configuration& config = mc->configuration();
    ref->set_model_params(config);
    for (int site_type = 0; site_type < config.num_site_types(); ++site_type) {
      ref->set_model_param("cutoff", site_type, cutoff_);
    }
    mc->add_to_reference(ref);
    if (use_cell_) {
      ref->set_visit_model_(MakeVisitModelCell({{"min_length",
                                                 str(cutoff_)}}));
    }
  } else {
    ASSERT(!use_cell_, "use_cell requires cutoff");
  }
}

}  // namespace feasst
