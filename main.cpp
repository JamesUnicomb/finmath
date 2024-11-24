#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <cassert>
#define assertm(exp, msg) assert(((void)msg, exp))

#include "tridiag.hpp"
#include "pde.hpp"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(optionpricing, mod)
{
    mod.doc() = "";

    mod.def("tridiag_solve", py::overload_cast<double, double, double, vector<double>&>(&tridiag_solve), "solve tridiagonal system.");
    mod.def("tridiag_mult", py::overload_cast<double, double, double, vector<double>&>(&tridiag_mult), "multiply a vector by a tridiagonal matrix.");

    py::class_<PDESolver>(mod, "PDESolver")
        .def(py::init<double, double, double>())
        .def("step", &PDESolver::step)
        .def_readwrite("alpha", &PDESolver::alpha)
        .def_readwrite("beta", &PDESolver::beta)
        .def_readwrite("gamma", &PDESolver::gamma);

#ifdef VERSION_INFO
    mod.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    mod.attr("__version__") = "dev";
#endif
}