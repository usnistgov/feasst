#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/utils.h"  // resize and fill
#include "utils/include/serialize_extra.h"
#include "utils/include/io.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "math/include/golden_search.h"
#include "math/include/formula.h"
#include "math/include/random_mt19937.h"
#include "math/include/utils_math.h"
#include "shape/include/shape.h"
#include "shape/include/shape_file.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select.h"
#include "confinement/include/model_table_cartesian_1d.h"

// HWH beware circular dependency
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

FEASST_MAPPER(ModelTableCart1D,);

ModelTableCart1D::ModelTableCart1D(argtype * args) {
  class_name_ = "ModelTableCart1D";
  dimension_ = integer("dimension", args, 2);
  table_file_ = str("table_file", args, "");
}
ModelTableCart1D::ModelTableCart1D(argtype args) : ModelTableCart1D(&args) {
  feasst_check_all_used(args);
}
ModelTableCart1D::~ModelTableCart1D() {}

void precompute_table1D(const Configuration& config, std::string table_file,
    std::vector<std::unique_ptr<Table1D> > * tables) {
  if (static_cast<int>(tables->size()) != 0) {
    return;
  }
  tables->resize(config.num_site_types());
  DEBUG("Read table");
  std::ifstream file(table_file);
  ASSERT(file.good(), "cannot find " << table_file);
  std::string line;
  std::getline(file, line);
  const std::vector<std::string> types = split(line, '=');
  ASSERT(types[0] == "site_types", "format error: " << types[0]);
  ASSERT(static_cast<int>(types.size()) == 2, "Error in formating of file:" <<
    table_file);
  for (const std::string& type : split(types[1], ',')) {
    const int index = config.site_type_name_to_index(type);
    std::getline(file, line);
    std::vector<std::string> values = split(line, ' ');
    const int num_values = values.size();
    (*tables)[index] = std::make_unique<Table1D>(argtype({{"num", str(num_values)}}));
    for (int v = 0; v < num_values; ++v) {
      (*tables)[index]->set_data(v, str_to_double(values[v]));
    }
  }
  std::stringstream errmsg;
  errmsg << "Extra lines in file:" << table_file << " . The number of lines "
    << "should be the number of site_types + 1.";
  if (!file.eof()) {
    std::getline(file, line);
    ASSERT(line.empty(), errmsg.str());
  }
  ASSERT(file.eof(), errmsg.str());
}

double ModelTableCart1D::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  TRACE(tables_.size());
  TRACE(site.type());
  const Table1D* table = tables_[site.type()].get();
  if (!table) {
    return 0.;
  }
  const double length = config.domain().side_length(dimension_);
  double val = wrapped_site.coord(dimension_)/length + 0.5;
  return table->forward_difference_interpolation(val);
}

void ModelTableCart1D::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(7890, ostr);
  feasst_serialize(dimension_, ostr);
  feasst_serialize(tables_, ostr);
}

ModelTableCart1D::ModelTableCart1D(std::istream& istr) : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 7890 && version <= 7890, "unrecognized verison: " << version);
  feasst_deserialize(&dimension_, istr);
  feasst_deserialize(&tables_, istr);
}

}  // namespace feasst
