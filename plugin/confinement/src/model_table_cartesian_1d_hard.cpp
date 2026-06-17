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
#include "confinement/include/model_table_cartesian_1d_hard.h"

// HWH beware circular dependency
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

FEASST_MAPPER(ModelTableCart1DHard, MakeTable1D());

ModelTableCart1DHard::ModelTableCart1DHard(std::shared_ptr<Table1D> table) {
  class_name_ = "ModelTableCart1DHard";
  tables_.resize(1);
  tables_[0] = table;
}

ModelTableCart1DHard::ModelTableCart1DHard(
    std::vector<std::shared_ptr<Table1D> > tables) {
  class_name_ = "ModelTableCart1DHard";
  tables_ = tables;
}

void ModelTableCart1DHard::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(9426, ostr);
  feasst_serialize(tables_, ostr);
}

ModelTableCart1DHard::ModelTableCart1DHard(std::istream& istr)
  : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9426, "unrecognized verison: " << version);
  int dim1;
  istr >> dim1;
  tables_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    // HWH for unk..
    //feasst_deserialize(table_, istr);
    int existing;
    istr >> existing;
    if (existing != 0) {
      tables_[index] = std::make_shared<Table1D>(istr);
    }
  }
}

double ModelTableCart1DHard::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  //const int type = site.type();
  //const double sigma = model_params.sigma().value(type);
  const std::vector<double>& sides = config.domain().side_lengths().coord();
  double val0 = 2.*wrapped_site.coord(0)/sides[0];
  if (val0 < 0) val0 *= -1;
  double val1 = 2.*wrapped_site.coord(1)/sides[1];
  if (val1 < 0) val1 *= -1;
  const double val1_max = tables_[site.type()]->linear_interpolation(val0);
  if (val1 > val1_max) return NEAR_INFINITY;
  return 0.;
}

const Table1D& ModelTableCart1DHard::table(const int site_type) const {
  return const_cast<Table1D&>(*tables_[site_type]);
}

class Objective1DHard : public Formula {
 public:
  Objective1DHard(Shape * shape, Position * point, const double diameter) {
    point_ = point, shape_ = shape; diameter_ = diameter;
  }
  double evaluate(const double x) const override {
    point_->set_coord(1, x);
    if (shape_->is_inside(*point_, diameter_)) return x;
    const double d = shape_->nearest_distance(*point_);
    return std::pow(d - diameter_, 2);
  }
 private:
  Shape * shape_;
  Position * point_;
  double diameter_;
};

void ModelTableCart1DHard::compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    argtype args,
    const int site_type) {
  const double diameter = dble("diameter", &args, 1.);
  Table1D * table = tables_[site_type].get();
  auto report = MakeProgressReport({{"num", str(table->num())}});
  Position point(domain->dimension());
  GoldenSearch minimize({{"tolerance", str(1e-8)},
    {"lower", "0"},
    {"upper", str(domain->side_length(0)/2)}});
  for (int bin = 0; bin < table->num(); ++bin) {
    point.set_coord(0, table->bin_value(bin)*domain->side_length(0)/2);
    Objective1DHard objective(shape, &point, diameter);
    table->set_data(bin, minimize.minimum(&objective));
    report->check();
  }
  feasst_check_all_used(args);
}

}  // namespace feasst
