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
#include "confinement/include/model_table_cartesian_2d_integr.h"

// HWH beware circular dependency
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

FEASST_MAPPER(ModelTableCart2DIntegr, MakeTable2D());

ModelTableCart2DIntegr::ModelTableCart2DIntegr(std::shared_ptr<Table2D> table) {
  class_name_ = "ModelTableCart2DIntegr";
  tables_.resize(1);
  tables_[0] = table;
}

ModelTableCart2DIntegr::ModelTableCart2DIntegr(
    std::vector<std::shared_ptr<Table2D> > tables) {
  class_name_ = "ModelTableCart2DIntegr";
  tables_ = tables;
}

void ModelTableCart2DIntegr::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(8460, ostr);
  feasst_serialize(tables_, ostr);
}

ModelTableCart2DIntegr::ModelTableCart2DIntegr(std::istream& istr)
  : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8460, "unrecognized verison: " << version);
  int dim1;
  istr >> dim1;
  tables_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    // HWH for unk..
    //feasst_deserialize(table_, istr);
    int existing;
    istr >> existing;
    if (existing != 0) {
      tables_[index] = std::make_shared<Table2D>(istr);
    }
  }
}

double ModelTableCart2DIntegr::energy(
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
  return tables_[site.type()]->linear_interpolation(val0, val1);
}

const Table2D& ModelTableCart2DIntegr::table(const int site_type) const {
  return const_cast<Table2D&>(*tables_[site_type]);
}

void ModelTableCart2DIntegr::compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    argtype integration_args,
    const int site_type) {
  Table2D * table = tables_[site_type].get();
  auto report = MakeProgressReport(
    {{"num", str(table->num0()*table->num1())}});
  Position point(domain->dimension());
  for (int bin0 = 0; bin0 < table->num0(); ++bin0) {
    point.set_coord(0,
      table->bin_to_value(0, bin0)*domain->side_length(0)/2);
    for (int bin1 = 0; bin1 < table->num1(); ++bin1) {
      point.set_coord(1,
        table->bin_to_value(1, bin1)*domain->side_length(1)/2.);
      if (shape->is_inside(point)) {
        const double en = -1*shape->integrate(point, random, integration_args);
        table->set_data(bin0, bin1, minimum(en, std::numeric_limits<float>::max()/1e10));
      }
      report->check();
    }
  }
}

void ModelTableCart2DIntegr::compute_table_omp(
    Shape * shape,
    Domain * domain,
    Random * random,
    argtype integration_args,
    const int site_type,
    const int node,
    const int num_nodes) {
  // allow shape internals to cache before parallel.
  shape->integrate(*MakePosition({{0., 0., 0.}}), random, integration_args);
  #ifdef _OPENMP
  Table2D * table = tables_[site_type].get();
  #pragma omp parallel
  {
    std::shared_ptr<Shape> shape_t = deep_copy_derived(shape);
    Domain domain_t = deep_copy(*domain);
    std::shared_ptr<Random> random_t = deep_copy_derived(random);
    auto thread = MakeThreadOMP();
    int iteration = 0;
    const int total = table->num0()*table->num1();
    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
    Position point(domain_t.dimension());
    for (int bin0 = 0; bin0 < table->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table->num1(); ++bin1) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        point.set_coord(0,
          table->bin_to_value(0, bin0)*domain_t.side_length(0)/2.);
        point.set_coord(1,
          table->bin_to_value(1, bin1)*domain_t.side_length(1)/2.);
        if (shape_t->is_inside(point)) {
          const double en = shape_t->integrate(point, random_t.get(), integration_args);
          table->set_data(bin0, bin1, minimum(en, std::numeric_limits<float>::max()/1e10));
        }
        report->check();
      }
      ++iteration;
    }}
  }
  #else // _OPENMP
    WARN("OMP not detected");
    compute_table(shape, domain, random, integration_args);
  #endif // _OPENMP
}

}  // namespace feasst
