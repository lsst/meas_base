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

#include "lsst/meas/base/SdssCentroid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyCentroidAlgorithm =
        py::class_<SdssCentroidAlgorithm, std::shared_ptr<SdssCentroidAlgorithm>, SimpleAlgorithm>;
using PyCentroidControl = py::class_<SdssCentroidControl>;
using PyCentroidTransform =
        py::class_<SdssCentroidTransform, std::shared_ptr<SdssCentroidTransform>, BaseTransform>;

PyCentroidControl declareCentroidControl(py::module &mod) {
    PyCentroidControl cls(mod, "SdssCentroidControl");

    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, binmax);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, peakMin);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, wfac);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, doFootprintCheck);
    LSST_DECLARE_CONTROL_FIELD(cls, SdssCentroidControl, maxDistToPeak);

    cls.def(py::init<>());

    return cls;
}

PyCentroidAlgorithm declareCentroidAlgorithm(py::module &mod) {
    PyCentroidAlgorithm cls(mod, "SdssCentroidAlgorithm");

    cls.attr("FAILURE") = py::cast(SdssCentroidAlgorithm::FAILURE);
    cls.attr("EDGE") = py::cast(SdssCentroidAlgorithm::EDGE);
    cls.attr("NO_SECOND_DERIVATIVE") = py::cast(SdssCentroidAlgorithm::NO_SECOND_DERIVATIVE);
    cls.attr("ALMOST_NO_SECOND_DERIVATIVE") = py::cast(SdssCentroidAlgorithm::ALMOST_NO_SECOND_DERIVATIVE);
    cls.attr("NOT_AT_MAXIMUM") = py::cast(SdssCentroidAlgorithm::NOT_AT_MAXIMUM);

    cls.def(py::init<SdssCentroidAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("measure", &SdssCentroidAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &SdssCentroidAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);

    return cls;
}

PyCentroidTransform declareCentroidTransform(py::module &mod) {
    PyCentroidTransform cls(mod, "SdssCentroidTransform");

    cls.def(py::init<SdssCentroidTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

}  // namespace

PYBIND11_MODULE(sdssCentroid, mod) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.transform");

    auto clsCentroidControl = declareCentroidControl(mod);
    auto clsCentroidAlgorithm = declareCentroidAlgorithm(mod);
    auto clsCentroidTransform = declareCentroidTransform(mod);

    clsCentroidAlgorithm.attr("Control") = clsCentroidControl;
    clsCentroidTransform.attr("Control") = clsCentroidControl;

    python::declareAlgorithm<SdssCentroidAlgorithm, SdssCentroidControl, SdssCentroidTransform>(
            clsCentroidAlgorithm, clsCentroidControl, clsCentroidTransform);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
