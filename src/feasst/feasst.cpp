#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
//#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/arguments_extra.h"
#include "utils/include/io.h"
#include "math/include/accumulator.h"
#include "math/include/position.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/log.h"
#include "actions/include/run.h"
#include "shape/include/shape.h"
#include "shape/include/half_space.h"
#include "shape/include/slab_sine.h"
#include "confinement/include/model_table_cartesian_3d_integr.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

void parse(feasst::MonteCarlo * mc, const std::string& line, const bool silent) {
  //std::cout << "py:feasst.cpp:line: " << line << std::endl;
  std::vector<std::string> lines = feasst::split(line, '\n');
  for (const std::string& line : lines) {
    if (!line.empty() && line[0] != '#') {
      auto parsed = feasst::parse_line(line, NULL, NULL);
      feasst::arglist list;
      list.push_back(parsed);
      mc->begin(list, silent);
    }
  }
}

//feasst::MonteCarlo deserialize(std::string str) {
//  std::stringstream ss;
//  ss << str;
//  return std::move(feasst::MonteCarlo(ss));
//}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 feasst plugin
        -----------------------

        .. currentmodule:: feasst

        .. autosummary::
           :toctree: _generate
    )pbdoc";
    m.def("parse", &parse,
        py::arg("mc"),
        py::arg("line"),
        py::arg("silent") = false,
        py::call_guard<py::gil_scoped_release>()); // required for concurrent threading

    py::class_<feasst::Accumulator>(m, "Accumulator")
        .def(py::init<>())
        .def(py::init<feasst::argtype>())
        .def("accumulate", &feasst::Accumulator::accumulate)
        .def("max_block_operations", &feasst::Accumulator::max_block_operations)
        .def("str", &feasst::Accumulator::str)
        .def("average", &feasst::Accumulator::average)
        .def("block_stdev", py::overload_cast<int, int>(&feasst::Accumulator::block_stdev, py::const_));

    py::class_<feasst::Position>(m, "Position")
        .def(py::init<feasst::argtype>())
        .def("coord", py::overload_cast<int>(&feasst::Position::coord, py::const_));

    py::class_<feasst::Random> random(m, "Random");
    random.def("uniform", static_cast<double (feasst::Random::*)()>(&feasst::Random::uniform));

    py::class_<feasst::RandomMT19937>(m, "RandomMT19937", random)
        .def(py::init<>())
        .def(py::init<feasst::argtype>());

    py::class_<feasst::MonteCarlo>(m, "MonteCarlo")
        .def(py::init<>())
        .def("configuration", &feasst::MonteCarlo::configuration)
        .def("is_criteria_set", &feasst::MonteCarlo::is_criteria_set)
        .def("serialize", py::overload_cast<>(&feasst::MonteCarlo::serialize, py::const_))
        .def("deserialize", &feasst::MonteCarlo::deserialize)
        ;

    py::class_<feasst::Configuration>(m, "Configuration")
        .def("particle", py::overload_cast<int>(&feasst::Configuration::particle, py::const_));

    py::class_<feasst::Particle>(m, "Particle")
        .def("site", &feasst::Particle::site);

    py::class_<feasst::Site>(m, "Site")
        .def("position", py::overload_cast<int>(&feasst::Site::position, py::const_));

    m.def("parse", &parse, R"pbdoc(
        Parse a single line with the text interface format.
    )pbdoc");

    py::class_<feasst::Shape> shape(m, "Shape");
    shape.def("integrate", static_cast<double (feasst::Shape::*)(const feasst::Position&, feasst::Random*, feasst::argtype)>(&feasst::Shape::integrate));

    py::class_<feasst::ShapeIntersect> shape_intersect(m, "ShapeIntersect", shape);

    py::class_<feasst::SlabSine>(m, "SlabSine", shape_intersect)
        .def(py::init<>())
        .def(py::init<feasst::argtype>());

    py::class_<feasst::HalfSpace>(m, "HalfSpace", shape)
        .def(py::init<>())
        .def(py::init<feasst::argtype>());

    py::class_<feasst::Domain>(m, "Domain")
        .def(py::init<feasst::argtype>())
        .def("side_length", &feasst::Domain::side_length);

    py::class_<feasst::Table3D>(m, "Table3D")
        .def(py::init<>())
        .def(py::init<feasst::argtype>())
        .def("data", &feasst::Table3D::data)
        .def("num", &feasst::Table3D::num)
        .def("minimum", &feasst::Table3D::minimum)
        .def("maximum", &feasst::Table3D::maximum)
        .def("value_to_nearest_bin", &feasst::Table3D::value_to_nearest_bin)
        .def("bin_to_value", &feasst::Table3D::bin_to_value)
        .def("serialize", py::overload_cast<>(&feasst::Table3D::serialize, py::const_))
        //.def("serialize", py::overload_cast<(feasst::Table3D::*)()>(&feasst::Table3D::serialize), py::const_)
        //.def("serialize", static_cast<std::string (feasst::Table3D::*)()>(&feasst::Table3D::serialize), py::const_)
        //.def("position", py::overload_cast<int>(&feasst::Site::position, py::const_));
        .def("deserialize", &feasst::Table3D::deserialize)
        ;

    py::class_<feasst::ModelTableCart3DIntegr>(m, "ModelTableCart3DIntegr")
        .def(py::init<feasst::argtype>())
        .def(py::init<std::shared_ptr<feasst::Table3D> >())
        .def("table", &feasst::ModelTableCart3DIntegr::table)
        .def("write", &feasst::ModelTableCart3DIntegr::write)
        .def("compute_table_omp", static_cast<void (feasst::ModelTableCart3DIntegr::*)(feasst::Shape*, const feasst::Domain&, feasst::Random*, feasst::argtype, const int, const int, const int)>(&feasst::ModelTableCart3DIntegr::compute_table_omp))
        .def("table", &feasst::ModelTableCart3DIntegr::table);

    py::class_<feasst::CheckEnergy>(m, "CheckEnergy")
        .def(py::init<feasst::argtype>());

    py::class_<feasst::Log>(m, "Log")
        .def(py::init<>())
        .def(py::init<feasst::argtype>())
        .def("write", static_cast<std::string (feasst::Log::*)(const feasst::MonteCarlo&)>(&feasst::Log::write))
        ;

    py::class_<feasst::Run>(m, "Run")
        .def(py::init<feasst::argtype>());

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
