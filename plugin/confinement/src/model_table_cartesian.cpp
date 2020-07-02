#include "utils/include/utils.h"  // resize and fill
#include "utils/include/serialize.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "math/include/golden_search.h"
#include "math/include/formula.h"
#include "shape/include/shape.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "confinement/include/model_table_cartesian.h"

// HWH beware circular dependency
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

class MapModelTableCart1DHard {
 public:
  MapModelTableCart1DHard() {
    auto model = MakeModelTableCart1DHard(MakeTable1D());
    model->deserialize_map()["ModelTableCart1DHard"] = model;
  }
};

static MapModelTableCart1DHard map_model_table_cart1d_ = MapModelTableCart1DHard();

void ModelTableCart1DHard::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(9426, ostr);
  feasst_serialize(table_, ostr);
}

ModelTableCart1DHard::ModelTableCart1DHard(std::istream& istr)
  : ModelOneBody() {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9426, "unrecognized verison: " << version);
  // HWH for unk..
  //feasst_deserialize(table_, istr);
  int existing;
  istr >> existing;
  if (existing != 0) {
    table_ = std::make_shared<Table1D>(istr);
  }
}

double ModelTableCart1DHard::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  //const int type = site.type();
  //const double sigma = model_params.sigma().value(type);
  const std::vector<double>& sides = config.domain().side_lengths().coord();
  double val0 = 2.*wrapped_site.coord(0)/sides[0];
  if (val0 < 0) val0 *= -1;
  double val1 = 2.*wrapped_site.coord(1)/sides[1];
  if (val1 < 0) val1 *= -1;
  const double val1_max = table_->linear_interpolation(val0);
  if (val1 > val1_max) return NEAR_INFINITY;
  return 0.;
}

const Table1D& ModelTableCart1DHard::table() const {
  return const_cast<Table1D&>(*table_);
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
    const argtype& args) {
  Arguments args_(args);
  const double diameter = args_.key("diameter").dflt("1.").dble();
  auto report = MakeProgressReport({{"num", str(table_->num())}});
  Position point(domain->dimension());
  GoldenSearch minimize({{"tolerance", str(1e-8)},
    {"lower", "0"},
    {"upper", str(domain->side_length(0)/2)}});
  for (int bin = 0; bin < table_->num(); ++bin) {
    point.set_coord(0, table_->bin_to_value(bin)*domain->side_length(0)/2);
    Objective1DHard objective(shape, &point, diameter);
    table_->set_data(bin, minimize.minimum(&objective));
    report->check();
  }
}

//void ModelTableCart1DHard::compute_table_omp(
//    Shape * shape,
//    Domain * domain,
//    Random * random,
//    const argtype& integration_args,
//    const int node,
//    const int num_nodes) {
//  // allow shape internals to cache before parallel.
//  shape->integrate(*MakePosition({{0., 0., 0.}}), random, integration_args);
//  #ifdef _OPENMP
//  #pragma omp parallel
//  {
//    std::shared_ptr<Shape> shape_t = deep_copy_derived(shape);
//    Domain domain_t = deep_copy(*domain);
//    std::shared_ptr<Random> random_t = deep_copy_derived(random);
//    auto thread = MakeThreadOMP();
//    int iteration = 0;
//    const int total = table_->num();
//    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
//    Position point(domain_t.dimension());
//    for (int bin0 = 0; bin0 < table_->num(); ++bin0) {
//      if (thread->in_chunk(iteration, total, node, num_nodes)) {
//        point.set_coord(0,
//          table_->bin_to_value(bin0)*domain_t.side_length(0)/2);
//        if (shape_t->is_inside(point)) {
//          table_->set_data(bin0,
//            shape_t->integrate(point, random_t.get(), integration_args));
//        }
//        report->check();
//      }
//      ++iteration;
//    }
//  }
//  #else // _OPENMP
//    WARN("OMP not detected");
//    compute_table(shape, domain, random, integration_args);
//  #endif // _OPENMP
//}

class MapModelTableCart2DIntegr {
 public:
  MapModelTableCart2DIntegr() {
    auto model = MakeModelTableCart2DIntegr(MakeTable2D());
    model->deserialize_map()["ModelTableCart2DIntegr"] = model;
  }
};

static MapModelTableCart2DIntegr map_model_table_cart2d_ = MapModelTableCart2DIntegr();

void ModelTableCart2DIntegr::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(8460, ostr);
  feasst_serialize(table_, ostr);
}

ModelTableCart2DIntegr::ModelTableCart2DIntegr(std::istream& istr)
  : ModelOneBody() {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8460, "unrecognized verison: " << version);
  // HWH for unk..
  //feasst_deserialize(table_, istr);
  int existing;
  istr >> existing;
  if (existing != 0) {
    table_ = std::make_shared<Table2D>(istr);
  }
}

double ModelTableCart2DIntegr::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  //const int type = site.type();
  //const double sigma = model_params.sigma().value(type);
  const std::vector<double>& sides = config.domain().side_lengths().coord();
  double val0 = 2.*wrapped_site.coord(0)/sides[0];
  if (val0 < 0) val0 *= -1;
  double val1 = 2.*wrapped_site.coord(1)/sides[1];
  if (val1 < 0) val1 *= -1;
  return table_->linear_interpolation(val0, val1);
}

const Table2D& ModelTableCart2DIntegr::table() const {
  return const_cast<Table2D&>(*table_);
}

void ModelTableCart2DIntegr::compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    const argtype& integration_args) {
  auto report = MakeProgressReport(
    {{"num", str(table_->num0()*table_->num1())}});
  Position point(domain->dimension());
  for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    point.set_coord(0,
      table_->bin_to_value(0, bin0)*domain->side_length(0)/2);
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
      point.set_coord(1,
        table_->bin_to_value(1, bin1)*domain->side_length(1)/2.);
      if (shape->is_inside(point)) {
        table_->set_data(bin0, bin1,
          -1*shape->integrate(point, random, integration_args));
      }
      report->check();
    }
  }
}

void ModelTableCart2DIntegr::compute_table_omp(
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
    const int total = table_->num0()*table_->num1();
    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
    Position point(domain_t.dimension());
    for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        point.set_coord(0,
          table_->bin_to_value(0, bin0)*domain_t.side_length(0)/2.);
        point.set_coord(1,
          table_->bin_to_value(1, bin1)*domain_t.side_length(1)/2.);
        if (shape_t->is_inside(point)) {
          table_->set_data(bin0, bin1,
            shape_t->integrate(point, random_t.get(), integration_args));
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

class MapModelTableCart3DIntegr {
 public:
  MapModelTableCart3DIntegr() {
    auto model = MakeModelTableCart3DIntegr(MakeTable3D());
    model->deserialize_map()["ModelTableCart3DIntegr"] = model;
  }
};

static MapModelTableCart3DIntegr map_model_table_cart3d_ = MapModelTableCart3DIntegr();

void ModelTableCart3DIntegr::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(6937, ostr);
  feasst_serialize(table_, ostr);
}

ModelTableCart3DIntegr::ModelTableCart3DIntegr(std::istream& istr)
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

double ModelTableCart3DIntegr::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  //const int type = site.type();
  //const double sigma = model_params.sigma().value(type);
  const std::vector<double>& sides = config.domain().side_lengths().coord();
  double val0 = 2.*wrapped_site.coord(0)/sides[0];
  if (val0 < 0) val0 *= -1;
  double val1 = 2.*wrapped_site.coord(1)/sides[1];
  if (val1 < 0) val1 *= -1;
  double val2 = 2.*wrapped_site.coord(2)/sides[2];
  if (val2 < 0) val2 *= -1;
  const double energy = table_->linear_interpolation(val0, val1, val2);
//  if (energy < -12) {
//    INFO(energy);
//    INFO(wrapped_site.str());
//    FATAL("err");
//  }
  return energy;
}

const Table3D& ModelTableCart3DIntegr::table() const {
  return const_cast<Table3D&>(*table_);
}

void ModelTableCart3DIntegr::compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    const argtype& integration_args) {
  auto report = MakeProgressReport(
    {{"num", str(table_->num0()*table_->num1()*table_->num2())}});
  Position point(domain->dimension());
  for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    point.set_coord(0,
      table_->bin_to_value(0, bin0)*domain->side_length(0)/2);
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
      point.set_coord(1,
        table_->bin_to_value(1, bin1)*domain->side_length(1)/2.);
      for (int bin2 = 0; bin2 < table_->num2(); ++bin2) {
        point.set_coord(2,
          table_->bin_to_value(2, bin2)*domain->side_length(2)/2.);
        if (shape->is_inside(point)) {
          table_->set_data(bin0, bin1, bin2,
            shape->integrate(point, random, integration_args));
        }
        report->check();
      }
    }
  }
}

void ModelTableCart3DIntegr::compute_table_omp(
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
    Position point(domain_t.dimension());
    for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
    for (int bin2 = 0; bin2 < table_->num2(); ++bin2) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        point.set_coord(0,
          table_->bin_to_value(0, bin0)*domain_t.side_length(0)/2);
        point.set_coord(1,
          table_->bin_to_value(1, bin1)*domain_t.side_length(1)/2.);
        point.set_coord(2,
          table_->bin_to_value(2, bin2)*domain_t.side_length(2)/2.);
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

void ModelTableCart3DIntegr::compute_table(
    System * system,
    Select * select) {
  auto report = MakeProgressReport(
    {{"num", str(table_->num0()*table_->num1()*table_->num2())}});
  ASSERT(select->num_sites() == 1, "assumes a single site");

  // MC machinery
  TrialSelect tsel;
  tsel.set_mobile(*select);
  PerturbAnywhere perturb;
  perturb.precompute(&tsel, system);

  const Domain& domain = system->configuration().domain();
  Position point(domain.dimension());
  for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    point.set_coord(0,
      table_->bin_to_value(0, bin0)*domain.side_length(0)/2);
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
      point.set_coord(1,
        table_->bin_to_value(1, bin1)*domain.side_length(1)/2.);
      for (int bin2 = 0; bin2 < table_->num2(); ++bin2) {
        point.set_coord(2,
          table_->bin_to_value(2, bin2)*domain.side_length(2)/2.);
        perturb.set_position(point, system, &tsel);
        system->energy();
        const double energy = system->perturbed_energy(*select);
        TRACE(system->configuration().select_particle(0).site(0).position().str() << " " << energy);
        table_->set_data(bin0, bin1, bin2, energy);
        perturb.finalize(system);
        system->finalize(*select);
        report->check();
      }
    }
  }
}

void ModelTableCart3DIntegr::compute_table_omp(
    System * system,
    Select * select,
    const int node,
    const int num_nodes) {
  #ifdef _OPENMP
  #pragma omp parallel
  {
    System system_t = deep_copy(*system);
    Select select_t = deep_copy(*select);

    // MC machinery
    TrialSelect tsel;
    tsel.set_mobile(select_t);
    PerturbAnywhere perturb;
    perturb.precompute(&tsel, &system_t);

    auto thread = MakeThreadOMP();
    int iteration = 0;
    const int total = table_->num0()*table_->num1()*table_->num2();
    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
    const Domain& domain = system->configuration().domain();
    Position point(domain.dimension());
    for (int bin0 = 0; bin0 < table_->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table_->num1(); ++bin1) {
    for (int bin2 = 0; bin2 < table_->num2(); ++bin2) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        point.set_coord(0,
          table_->bin_to_value(0, bin0)*domain.side_length(0)/2);
        point.set_coord(1,
          table_->bin_to_value(1, bin1)*domain.side_length(1)/2.);
        point.set_coord(2,
          table_->bin_to_value(2, bin2)*domain.side_length(2)/2.);

        perturb.set_position(point, &system_t, &tsel);
        system_t.energy();
        const double energy = system_t.perturbed_energy(select_t);
        TRACE(system_t.configuration().select_particle(0).site(0).position().str() << " " << energy);
        table_->set_data(bin0, bin1, bin2, energy);
        perturb.finalize(&system_t);
        system_t.finalize(select_t);

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
