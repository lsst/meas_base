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
#include "pybind11/stl.h"

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/PsfFlux.h"
#include "lsst/meas/base/FluxUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm = py::class_<PsfFluxAlgorithm, std::shared_ptr<PsfFluxAlgorithm>, SimpleAlgorithm>;
using PyFluxControl = py::class_<PsfFluxControl>;
using PyFluxTransform = py::class_<PsfFluxTransform, std::shared_ptr<PsfFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "PsfFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, PsfFluxControl, badMaskPlanes);

    cls.def(py::init<>());

    return cls;
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "PsfFluxAlgorithm");

    cls.attr("FAILURE") = py::cast(PsfFluxAlgorithm::FAILURE);
    cls.attr("NO_GOOD_PIXELS") = py::cast(PsfFluxAlgorithm::NO_GOOD_PIXELS);
    cls.attr("EDGE") = py::cast(PsfFluxAlgorithm::EDGE);

    cls.def(py::init<PsfFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    return cls;
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "PsfFluxTransform");

    cls.def(py::init<PsfFluxTransform::Control const &, std::string const &, afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

}  // <anonymous>

PYBIND11_PLUGIN(psfFlux) {
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.fluxUtilities");
    py::module::import("lsst.meas.base.transform");

    py::module mod("psfFlux");

    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<PsfFluxAlgorithm, PsfFluxControl, PsfFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
