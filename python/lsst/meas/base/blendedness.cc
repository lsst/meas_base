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
#include "lsst/utils/python.h"

#include <memory>

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/Blendedness.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyBlendednessAlgorithm =
        py::class_<BlendednessAlgorithm, std::shared_ptr<BlendednessAlgorithm>, SimpleAlgorithm>;
using PyBlendednessControl = py::class_<BlendednessControl>;

PyBlendednessControl declareBlendednessControl(lsst::utils::python::WrapperCollection &wrappers) {
   return wrappers.wrapType(PyBlendednessControl(wrappers.module, "BlendednessControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, doOld);
        LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, doFlux);
        LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, doShape);
        LSST_DECLARE_CONTROL_FIELD(cls, BlendednessControl, nSigmaWeightMax);

        cls.def(py::init<>());
    });
}

PyBlendednessAlgorithm declareBlendednessAlgorithm(lsst::utils::python::WrapperCollection &wrappers) {
    return wrappers.wrapType(PyBlendednessAlgorithm(wrappers.module, "BlendednessAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<BlendednessAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.attr("FAILURE") = py::cast(BlendednessAlgorithm::FAILURE);
        cls.attr("NO_CENTROID") = py::cast(BlendednessAlgorithm::NO_CENTROID);
        cls.attr("NO_SHAPE") = py::cast(BlendednessAlgorithm::NO_SHAPE);

        cls.def_static("computeAbsExpectation", &BlendednessAlgorithm::computeAbsExpectation, "data"_a,
                       "variance"_a);
        cls.def_static("computeAbsBias", &BlendednessAlgorithm::computeAbsBias, "mu"_a, "variance"_a);
        cls.def("measureChildPixels", &BlendednessAlgorithm::measureChildPixels, "image"_a, "child"_a);
        cls.def("measureParentPixels", &BlendednessAlgorithm::measureParentPixels, "image"_a, "child"_a);
        cls.def("measure", &BlendednessAlgorithm::measure, "measRecord"_a, "exposure"_a);
        cls.def("fail", &BlendednessAlgorithm::measure, "measRecord"_a, "error"_a = nullptr);
    });
}

}  // namespace

void wrapBlendedness(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.addInheritanceDependency("lsst.afw.table");
    // Depends on algorithm
    // Depends on flagHandler

    auto clsBlendednessControl = declareBlendednessControl(wrappers);
    auto clsBlendednessAlgorithm = declareBlendednessAlgorithm(wrappers);

    clsBlendednessAlgorithm.attr("Control") = clsBlendednessControl;

    python::declareAlgorithm<BlendednessAlgorithm, BlendednessControl>(clsBlendednessAlgorithm,
                                                                       clsBlendednessControl);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
