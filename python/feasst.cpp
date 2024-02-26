#include <pybind11/pybind11.h>
#include <string>
#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

void parse(feasst::MonteCarlo * mc, const std::string& line) {
  auto parsed = feasst::parse_line(line, NULL, NULL);
  feasst::arglist list;
  list.push_back(parsed);
  mc->begin(list);
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 feasst plugin
        -----------------------

        .. currentmodule:: feasst

        .. autosummary::
           :toctree: _generate
    )pbdoc";

    py::class_<feasst::MonteCarlo>(m, "MonteCarlo")
        .def(py::init<>());

    m.def("parse", &parse, R"pbdoc(
        Parse a single line with the text interface format.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
