/*
 * LSST Data Management System
 * Copyright 2008-2017  AURA/LSST.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "pybind11/pybind11.h"

#include <memory>

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/GaussianFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm =
        py::class_<GaussianFluxAlgorithm, std::shared_ptr<GaussianFluxAlgorithm>, SimpleAlgorithm>;
using PyFluxControl = py::class_<GaussianFluxControl>;
using PyFluxTransform =
        py::class_<GaussianFluxTransform, std::shared_ptr<GaussianFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "GaussianFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, GaussianFluxControl, background);

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "GaussianFluxAlgorithm");

    cls.def(py::init<GaussianFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.attr("FAILURE") = py::cast(GaussianFluxAlgorithm::FAILURE);

    cls.def("measure", &GaussianFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &GaussianFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "GaussianFluxTransform");

    cls.def(py::init<GaussianFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

}  // namespace

PYBIND11_PLUGIN(gaussianFlux) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.transform");

    py::module mod("gaussianFlux");

    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<GaussianFluxAlgorithm, GaussianFluxControl, GaussianFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);

    return mod.ptr();
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
