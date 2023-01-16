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
#include "pybind11/stl.h"
#include "lsst/utils/python.h"

#include <memory>

#include "lsst/pex/config/python.h"

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/PixelFlags.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyAlgorithm =
        py::class_<PixelFlagsAlgorithm, std::shared_ptr<PixelFlagsAlgorithm>, SimpleAlgorithm>;
using PyControl = py::class_<PixelFlagsControl>;

void declareControl(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyControl(wrappers.module, "PixelFlagsControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, PixelFlagsControl, masksFpAnywhere);
        LSST_DECLARE_CONTROL_FIELD(cls, PixelFlagsControl, masksFpCenter);

        cls.def(py::init<>());
    });
}

void declareAlgorithm(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(PyAlgorithm(wrappers.module, "PixelFlagsAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(
                py::init<PixelFlagsAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);

        cls.def("measure", &PixelFlagsAlgorithm::measure, "measRecord"_a, "exposure"_a);
        cls.def("fail", &PixelFlagsAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    });
}

}  // namespace

void wrapPixelFlags(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.addSignatureDependency("lsst.afw.table");
    // Depends on algorithm

    declareControl(wrappers);
    declareAlgorithm(wrappers);

    // The original code did not have a python::declareAlgorithm call.  Was this a bug?
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
