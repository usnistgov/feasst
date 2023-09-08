#include <fstream>
#include "utils/include/utils.h"  // resize and fill
#include "utils/include/serialize.h"
#include "utils/include/progress_report.h"
#include "threads/include/thread_omp.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "math/include/golden_search.h"
#include "math/include/formula.h"
#include "math/include/random_mt19937.h"
#include "shape/include/shape.h"
#include "shape/include/shape_file.h"
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
    point.set_coord(0, table->bin_to_value(bin)*domain->side_length(0)/2);
    Objective1DHard objective(shape, &point, diameter);
    table->set_data(bin, minimize.minimum(&objective));
    report->check();
  }
  FEASST_CHECK_ALL_USED(args);
}

class MapModelTableCart2DIntegr {
 public:
  MapModelTableCart2DIntegr() {
    auto model = MakeModelTableCart2DIntegr(MakeTable2D());
    model->deserialize_map()["ModelTableCart2DIntegr"] = model;
  }
};

static MapModelTableCart2DIntegr map_model_table_cart2d_ = MapModelTableCart2DIntegr();

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
        table->set_data(bin0, bin1,
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
          table->set_data(bin0, bin1,
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

ModelTableCart3DIntegr::ModelTableCart3DIntegr(std::shared_ptr<Table3D> table) {
  class_name_ = "ModelTableCart3DIntegr";
  tables_.resize(1);
  tables_[0] = table;
}

ModelTableCart3DIntegr::ModelTableCart3DIntegr(
    std::vector<std::shared_ptr<Table3D> > tables) {
  class_name_ = "ModelTableCart3DIntegr";
  tables_ = tables;
}

ModelTableCart3DIntegr::ModelTableCart3DIntegr(argtype *args) {
  class_name_ = "ModelTableCart3DIntegr";
  const std::string file_name = str("file_name", args);
  // if file_name is the only argument, then read from file
  INFO(args->size());
  if (static_cast<int>(args->size()) == 0) {
    read(file_name);
    return;
  }
  // otherwise, generate the table and write
  auto random = RandomMT19937().factory(str("random", args, "RandomMT19937"), args);
  INFO(random->class_name());
  auto domain = std::make_shared<Domain>(args);
  ASSERT(domain->dimension() != 0, "Domain arguments are required.");
  INFO("args " << str(*args));
  auto shape = MakeShapeFile({{"file_name", str("shape_file_name", args)}});
  INFO(shape->class_name());
  //ASSERT(config.num_site_types() <= 1, "only implemented for single site");
  const int site_type = 0;
  tables_.push_back(std::make_shared<Table3D>(args));
  site_types_.push_back(site_type);
  const bool use_omp = boolean("use_omp", args, false);
  INFO("use_omp " << use_omp);
  if (use_omp) {
    const int node = integer("node", args, 0);
    const int num_node = integer("num_node", args, 1);
    compute_table_omp(shape.get(), *domain, random.get(), args, site_type, node, num_node);
  } else {
    compute_table(shape.get(), *domain, random.get(), args, site_type);
  }
  write(file_name);
}
ModelTableCart3DIntegr::ModelTableCart3DIntegr(argtype args) : ModelTableCart3DIntegr(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void ModelTableCart3DIntegr::read(const std::string file_name) {
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find " << file_name);
  std::string descript, line;
  int int_val;
  file >> descript >> int_val;
  ASSERT(descript == "site_types", "format error: " << descript);
  const int num_sites = int_val;
  DEBUG("num_sites " << num_sites);
  site_types_.resize(num_sites);
  for (int type = 0; type < num_sites; ++type) {
    file >> int_val;
    DEBUG("site " << int_val);
    site_types_[type] = int_val;
  }
  std::getline(file, line);
  Table3D empty;
  for (int type = 0; type < num_sites; ++type) {
    std::getline(file, line);
    tables_.push_back(std::make_shared<Table3D>(empty.deserialize(line)));
  }
  tables_.resize(num_sites);
}

void ModelTableCart3DIntegr::write(const std::string file_name) const {
  std::ofstream file(file_name);
  const int num_site_types = static_cast<int>(site_types_.size());
  if (num_site_types != 0) {
    file << "site_types " << num_site_types << " ";
    for (int site : site_types_) {
      file << site << " ";
    }
  } else {
    file << "site_types 1 0";
  }
  file << std::endl;
  for (std::shared_ptr<Table3D> table : tables_) {
    file << table->serialize() << std::endl;
  }
}

void ModelTableCart3DIntegr::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(6937, ostr);
  feasst_serialize(tables_, ostr);
  feasst_serialize(site_types_, ostr);
}

ModelTableCart3DIntegr::ModelTableCart3DIntegr(std::istream& istr)
  : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6937, "unrecognized verison: " << version);
  int dim1;
  istr >> dim1;
  tables_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    // HWH for unk..
    //feasst_deserialize(table_, istr);
    int existing;
    istr >> existing;
    if (existing != 0) {
      tables_[index] = std::make_shared<Table3D>(istr);
    }
  }
  feasst_deserialize(&site_types_, istr);
}

double ModelTableCart3DIntegr::energy(
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
  double val2 = 2.*wrapped_site.coord(2)/sides[2];
  if (val2 < 0) val2 *= -1;
  TRACE(tables_.size());
  TRACE(site.type());
  return tables_[site.type()]->linear_interpolation(val0, val1, val2);
}

const Table3D& ModelTableCart3DIntegr::table(const int site_type) const {
  return const_cast<Table3D&>(*tables_[site_type]);
}

void ModelTableCart3DIntegr::compute_table(
    Shape * shape,
    const Domain& domain,
    Random * random,
    argtype * integration_args,
    const int site_type) {
  argtype integration_args_copy = *integration_args;
  // allow shape internals to cache and parse arguments
  shape->integrate(*MakePosition({{0., 0., 0.}}), random, integration_args);
  Table3D * table = tables_[site_type].get();
  auto report = MakeProgressReport(
    {{"num", str(table->num0()*table->num1()*table->num2())}});
  Position point(domain.dimension());
  for (int bin0 = 0; bin0 < table->num0(); ++bin0) {
    point.set_coord(0,
      table->bin_to_value(0, bin0)*domain.side_length(0)/2);
    for (int bin1 = 0; bin1 < table->num1(); ++bin1) {
      point.set_coord(1,
        table->bin_to_value(1, bin1)*domain.side_length(1)/2.);
      for (int bin2 = 0; bin2 < table->num2(); ++bin2) {
        point.set_coord(2,
          table->bin_to_value(2, bin2)*domain.side_length(2)/2.);
        if (shape->is_inside(point)) {
          table->set_data(bin0, bin1, bin2,
            shape->integrate(point, random, integration_args_copy));
        }
        report->check();
      }
    }
  }
}
void ModelTableCart3DIntegr::compute_table(
    Shape * shape,
    const Domain& domain,
    Random * random,
    argtype integration_args,
    const int site_type) {
  compute_table(shape, domain, random, &integration_args, site_type);
}

void ModelTableCart3DIntegr::compute_table_omp(
    Shape * shape,
    const Domain& domain,
    Random * random,
    argtype * integration_args,
    const int site_type,
    const int node,
    const int num_nodes) {
  argtype integration_args_copy = *integration_args;
  // allow shape internals to cache before parallel.
  shape->integrate(*MakePosition({{0., 0., 0.}}), random, integration_args);
  #ifdef _OPENMP
  Table3D * table = tables_[site_type].get();
  #pragma omp parallel
  {
    std::shared_ptr<Shape> shape_t = deep_copy_derived(shape);
    Domain domain_t = domain;
    std::shared_ptr<Random> random_t = deep_copy_derived(random);
    auto thread = MakeThreadOMP();
    int iteration = 0;
    const int total = table->num0()*table->num1()*table->num2();
    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
    Position point(domain_t.dimension());
    for (int bin0 = 0; bin0 < table->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table->num1(); ++bin1) {
    for (int bin2 = 0; bin2 < table->num2(); ++bin2) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        point.set_coord(0,
          table->bin_to_value(0, bin0)*domain_t.side_length(0)/2);
        point.set_coord(1,
          table->bin_to_value(1, bin1)*domain_t.side_length(1)/2.);
        point.set_coord(2,
          table->bin_to_value(2, bin2)*domain_t.side_length(2)/2.);
        if (shape_t->is_inside(point)) {
          table->set_data(bin0, bin1, bin2,
            shape_t->integrate(point, random_t.get(), integration_args_copy));
        }
        report->check();
      }
      ++iteration;
    }}}
  }
  #else // _OPENMP
    WARN("OMP not detected");
    compute_table(shape, domain, random, &integration_args_copy);
  #endif // _OPENMP
}

void ModelTableCart3DIntegr::compute_table_omp(
    Shape * shape,
    const Domain& domain,
    Random * random,
    argtype integration_args,
    const int site_type,
    const int node,
    const int num_nodes) {
  compute_table_omp(shape, domain, random, &integration_args, site_type, node, num_nodes);
}

void ModelTableCart3DIntegr::compute_table(
    System * system,
    Select * select,
    const int site_type) {

  // for each unique site_type in select, ... (or, use particle_type as input instead of select+site_type.. add/del part physical part from system during making of table...

  Table3D * table = tables_[site_type].get();
  auto report = MakeProgressReport(
    {{"num", str(table->num0()*table->num1()*table->num2())}});
  ASSERT(select->num_sites() == 1, "assumes a single site");

  // MC machinery
  TrialSelect tsel;
  tsel.set_mobile(*select);
  PerturbAnywhere perturb;
  perturb.precompute(&tsel, system);

  const Domain& domain = system->configuration().domain();
  Position point(domain.dimension());
  for (int bin0 = 0; bin0 < table->num0(); ++bin0) {
    point.set_coord(0,
      table->bin_to_value(0, bin0)*domain.side_length(0)/2);
    for (int bin1 = 0; bin1 < table->num1(); ++bin1) {
      point.set_coord(1,
        table->bin_to_value(1, bin1)*domain.side_length(1)/2.);
      for (int bin2 = 0; bin2 < table->num2(); ++bin2) {
        point.set_coord(2,
          table->bin_to_value(2, bin2)*domain.side_length(2)/2.);
        perturb.set_position(point, system, &tsel);
        // HWH configuration_index_
        system->energy(0);
        // HWH configuration_index_
        const double energy = system->perturbed_energy(*select, 0);
        TRACE(system->configuration().select_particle(0).site(0).position().str() << " " << energy);
        table->set_data(bin0, bin1, bin2, energy);
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
    const int site_type,
    const int node,
    const int num_nodes) {
  #ifdef _OPENMP
  Table3D * table = tables_[site_type].get();
  #pragma omp parallel
  {
    System system_t = deep_copy(*system);
    Select select_t = deep_copy(*select);

    // MC machinery
    TrialSelect tsel;
    tsel.precompute(&system_t);
    tsel.set_mobile(select_t);
    PerturbAnywhere perturb;
    perturb.precompute(&tsel, &system_t);

    auto thread = MakeThreadOMP();
    int iteration = 0;
    const int total = table->num0()*table->num1()*table->num2();
    auto report = MakeProgressReport({{"num", str(total/thread->num()/num_nodes)}});
    const Domain& domain = system->configuration().domain();
    Position point(domain.dimension());
    for (int bin0 = 0; bin0 < table->num0(); ++bin0) {
    for (int bin1 = 0; bin1 < table->num1(); ++bin1) {
    for (int bin2 = 0; bin2 < table->num2(); ++bin2) {
      if (thread->in_chunk(iteration, total, node, num_nodes)) {
        point.set_coord(0,
          table->bin_to_value(0, bin0)*domain.side_length(0)/2);
        point.set_coord(1,
          table->bin_to_value(1, bin1)*domain.side_length(1)/2.);
        point.set_coord(2,
          table->bin_to_value(2, bin2)*domain.side_length(2)/2.);

        perturb.set_position(point, &system_t, &tsel);
        perturb.set_finalize_possible(true, &tsel);
        // HWH configuration_index_
        system_t.energy(0);
        // HWH configuration_index_
        const double energy = system_t.perturbed_energy(select_t, 0);
        TRACE(system_t.configuration().select_particle(0).site(0).position().str() << " " << energy);
        table->set_data(bin0, bin1, bin2, energy);
        perturb.finalize(&system_t);
        system_t.finalize(select_t);

        report->check();
      }
      ++iteration;
    }}}
  }
  #else // _OPENMP
    WARN("OMP not detected");
    compute_table(system, select, site_type);
  #endif // _OPENMP
}

}  // namespace feasst
