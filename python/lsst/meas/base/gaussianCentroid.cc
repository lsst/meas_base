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

#include "lsst/meas/base/GaussianCentroid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyCentroidAlgorithm =
        py::class_<GaussianCentroidAlgorithm, std::shared_ptr<GaussianCentroidAlgorithm>, SimpleAlgorithm>;
using PyCentroidControl = py::class_<GaussianCentroidControl>;
using PyCentroidTransform =
        py::class_<GaussianCentroidTransform, std::shared_ptr<GaussianCentroidTransform>, CentroidTransform>;

PyCentroidControl declareCentroidControl(py::module & mod) {
    PyCentroidControl cls(mod, "GaussianCentroidControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, GaussianCentroidControl, doFootprintCheck);
    LSST_DECLARE_CONTROL_FIELD(cls, GaussianCentroidControl, maxDistToPeak);

    return cls;
}

template <typename Pixel>
void declareFitCentroid(PyCentroidAlgorithm & cls) {
    cls.def_static("fitCentroid", &GaussianCentroidAlgorithm::fitCentroid<Pixel>, "im"_a, "x0"_a, "y0"_a);
}

PyCentroidAlgorithm declareCentroidAlgorithm(py::module & mod) {
    PyCentroidAlgorithm cls(mod, "GaussianCentroidAlgorithm");

    cls.attr("FAILURE") = py::cast(GaussianCentroidAlgorithm::FAILURE);
    cls.attr("NO_PEAK") = py::cast(GaussianCentroidAlgorithm::NO_PEAK);

    cls.def(py::init<GaussianCentroidAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    declareFitCentroid<float>(cls);
    declareFitCentroid<double>(cls);

    cls.def("measure", &GaussianCentroidAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &GaussianCentroidAlgorithm::fail, "measRecord"_a, "error"_a=nullptr);

    return cls;
}

}  // <anonymous>

PYBIND11_PLUGIN(gaussianCentroid) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.transform");

    py::module mod("gaussianCentroid");

    auto clsCentroidControl = declareCentroidControl(mod);
    auto clsCentroidAlgorithm = declareCentroidAlgorithm(mod);
    PyCentroidTransform clsCentroidTransform(mod, "GaussianCentroidTransform");

    clsCentroidAlgorithm.attr("Control") = clsCentroidControl;
    clsCentroidTransform.attr("Control") = clsCentroidControl;

    python::declareAlgorithm<GaussianCentroidAlgorithm, GaussianCentroidControl, GaussianCentroidTransform>(
            clsCentroidAlgorithm, clsCentroidControl, clsCentroidTransform);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
