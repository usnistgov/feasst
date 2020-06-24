#include "utils/include/utils.h"  // resize and fill
#include "utils/include/serialize.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "shape/include/shape.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "confinement/include/model_table_cartesian.h"

namespace feasst {

class MapModelTableCart3FoldSym {
 public:
  MapModelTableCart3FoldSym() {
    auto model = MakeModelTableCart3FoldSym(MakeTable3D());
    model->deserialize_map()["ModelTableCart3FoldSym"] = model;
  }
};

static MapModelTableCart3FoldSym map_model_table_cartesian_ = MapModelTableCart3FoldSym();

void ModelTableCart3FoldSym::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(6937, ostr);
  feasst_serialize(table_, ostr);
}

ModelTableCart3FoldSym::ModelTableCart3FoldSym(std::istream& istr)
  : ModelOneBody() {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6937, "unrecognized verison: " << version);
  // HWH for unk..
  //feasst_deserialize(table_, istr);
  int existing;
  istr >> existing;
  if (existing != 0) {
    table_ = std::make_shared<Table3D>(istr);
  }
}

double ModelTableCart3FoldSym::energy(
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  //const int type = site.type();
  //const double sigma = model_params.sigma().value(type);
  const std::vector<double>& sides = config.domain().side_lengths().coord();
  double val0 = 2.*site.position().coord(0)/sides[0];
  if (val0 < 0) val0 *= -1;
  double val1 = 2.*site.position().coord(1)/sides[1];
  if (val1 < 0) val1 *= -1;
  double val2 = 2.*site.position().coord(2)/sides[2];
  if (val2 < 0) val2 *= -1;
  return table_->linear_interpolation(val0, val1, val2);
}

const Table3D& ModelTableCart3FoldSym::table() const {
  return const_cast<Table3D&>(*table_);
}

void ModelTableCart3FoldSym::compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    const argtype& integration_args) {
  auto report = MakeProgressReport(
    {{"num", str(table_->num0()*table_->num1()*table_->num2())}});
  for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
  for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
  for (int bin2 = 0; bin2 < table_->num2(); ++bin2) {
    Position point({
      table_->bin_to_value(0, bin0)*domain->side_length(0)/2,
      table_->bin_to_value(1, bin1)*domain->side_length(1)/2.,
      table_->bin_to_value(2, bin2)*domain->side_length(2)/2.});
    point.add(domain->shift_opt(point));
    if (shape->is_inside(point)) {
      table_->set_data(bin0, bin1, bin2,
        -1*shape->integrate(point, random, integration_args));
    }
    report->check();
  }}}
}

void ModelTableCart3FoldSym::compute_table_omp(
    Shape * shape,
    Domain * domain,
    Random * random,
    const argtype& integration_args,
    const int node,
    const int num_nodes) {
  // allow shape internals to cache before parallel.
  shape->integrate(*MakePosition({{0., 0., 0.}}), random, integration_args);
  #ifdef _OPENMP
  #pragma omp parallel
  {
    std::shared_ptr<Shape> shape_t = deep_copy_derived(shape);
    Domain domain_t = deep_copy(*domain);
    std::shared_ptr<Random> random_t = deep_copy_derived(random);
    auto thread = MakeThreadOMP();
    int iteration = 0;
    const int total = table_->num0()*table_->num1()*table_->num2();
    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
    for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
    for (int bin2 = 0; bin2 < table_->num2(); ++bin2) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        Position point({
          table_->bin_to_value(0, bin0)*domain_t.side_length(0)/2,
          table_->bin_to_value(1, bin1)*domain_t.side_length(1)/2.,
          table_->bin_to_value(2, bin2)*domain_t.side_length(2)/2.});
        point.add(domain_t.shift_opt(point));
        if (shape_t->is_inside(point)) {
          table_->set_data(bin0, bin1, bin2,
            shape_t->integrate(point, random_t.get(), integration_args));
        }
        report->check();
      }
      ++iteration;
    }}}
  }
  #else // _OPENMP
    WARN("OMP not detected");
    compute_table(shape, domain, random, integration_args);
  #endif // _OPENMP
}

}  // namespace feasst
