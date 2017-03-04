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

#include "lsst/meas/base/NaiveCentroid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyCentroidAlgorithm =
        py::class_<NaiveCentroidAlgorithm, std::shared_ptr<NaiveCentroidAlgorithm>, SimpleAlgorithm>;
using PyCentroidControl = py::class_<NaiveCentroidControl>;
using PyCentroidTransform =
        py::class_<NaiveCentroidTransform, std::shared_ptr<NaiveCentroidTransform>, CentroidTransform>;

PyCentroidControl declareCentroidControl(py::module &mod) {
    PyCentroidControl cls(mod, "NaiveCentroidControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, NaiveCentroidControl, background);
    LSST_DECLARE_CONTROL_FIELD(cls, NaiveCentroidControl, doFootprintCheck);
    LSST_DECLARE_CONTROL_FIELD(cls, NaiveCentroidControl, maxDistToPeak);

    return cls;
}

PyCentroidAlgorithm declareCentroidAlgorithm(py::module &mod) {
    PyCentroidAlgorithm cls(mod, "NaiveCentroidAlgorithm");

    cls.attr("FAILURE") = py::cast(NaiveCentroidAlgorithm::FAILURE);
    cls.attr("NO_COUNTS") = py::cast(NaiveCentroidAlgorithm::NO_COUNTS);
    cls.attr("EDGE") = py::cast(NaiveCentroidAlgorithm::EDGE);

    cls.def(py::init<NaiveCentroidAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &NaiveCentroidAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &NaiveCentroidAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyCentroidTransform declareCentroidTransform(py::module &mod) {
    PyCentroidTransform cls(mod, "NaiveCentroidTransform");

    cls.def(py::init<NaiveCentroidTransform::Control const &, std::string const &, afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

}  // <anonymous>

PYBIND11_PLUGIN(naiveCentroid) {
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.transform");

    py::module mod("naiveCentroid");

    /* Module level */
    auto clsCentroidControl = declareCentroidControl(mod);
    auto clsCentroidAlgorithm = declareCentroidAlgorithm(mod);
    auto clsCentroidTransform = declareCentroidTransform(mod);

    clsCentroidAlgorithm.attr("Control") = clsCentroidControl;
    clsCentroidTransform.attr("Control") = clsCentroidControl;

    /* Members */
    python::declareAlgorithm<NaiveCentroidAlgorithm, NaiveCentroidControl, NaiveCentroidTransform>(
            clsCentroidAlgorithm, clsCentroidControl, clsCentroidTransform);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
