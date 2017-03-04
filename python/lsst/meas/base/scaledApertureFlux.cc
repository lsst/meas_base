/* 
 * LSST Data Management System
 * Copyright 2008-2016  AURA/LSST.
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
#include <memory>

#include "pybind11/pybind11.h"

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/ScaledApertureFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm = py::class_<ScaledApertureFluxAlgorithm, std::shared_ptr<ScaledApertureFluxAlgorithm>,
                                   SimpleAlgorithm>;
using PyFluxControl = py::class_<ScaledApertureFluxControl>;
using PyFluxTransform =
        py::class_<ScaledApertureFluxTransform, std::shared_ptr<ScaledApertureFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "ScaledApertureFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, ScaledApertureFluxControl, scale);
    LSST_DECLARE_CONTROL_FIELD(cls, ScaledApertureFluxControl, shiftKernel);

    cls.def(py::init<>());

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "ScaledApertureFluxAlgorithm");

    cls.def(py::init<ScaledApertureFluxAlgorithm::Control const &, std::string const &,
                     afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &ScaledApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &ScaledApertureFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "ScaledApertureFluxTransform");

    cls.def(py::init<ScaledApertureFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

}  // <anonymous>

PYBIND11_PLUGIN(scaledApertureFlux) {
    py::module mod("scaledApertureFlux");

    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<ScaledApertureFluxAlgorithm, ScaledApertureFluxControl,
                             ScaledApertureFluxTransform>(clsFluxAlgorithm, clsFluxControl, clsFluxTransform);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
