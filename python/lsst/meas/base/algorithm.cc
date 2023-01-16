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

//#include "lsst/meas/base/python.h"

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/Algorithm.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

void declareBaseAlgorithm(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(
            py::class_<BaseAlgorithm, std::shared_ptr<BaseAlgorithm>>(wrappers.module,
                                                                      "BaseAlgorithm"),
            [](auto &mod, auto &cls) {
                cls.def("fail", &BaseAlgorithm::fail, "measRecord"_a, "error"_a = NULL);
                cls.def("getLogName", &SimpleAlgorithm::getLogName);
    });
}

void declareSingleFrameAlgorithm(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(
            py::class_<SingleFrameAlgorithm, std::shared_ptr<SingleFrameAlgorithm>>(wrappers.module,
                                                                                    "SingleFrameAlgorithm"),
            [](auto &mod, auto &cls) {
                cls.def("measure", &SingleFrameAlgorithm::measure, "record"_a, "exposure"_a);
    });
}

void declareSimpleAlgorithm(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(
            py::class_<SimpleAlgorithm, std::shared_ptr<SimpleAlgorithm>>(wrappers.module,
                                                                                    "SimpleAlgorithm"),
            [](auto &mod, auto &cls) {
                cls.def("measureForced", &SimpleAlgorithm::measureForced, "measRecord"_a, "exposure"_a,
                           "refRecord"_a, "refWcs"_a);
    });
}

}  // namespace

void wrapAlgorithm(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.addInheritanceDependency("lsst.afw.image");
    wrappers.addInheritanceDependency("lsst.afw.table");

    declareBaseAlgorithm(wrappers);
    declareSingleFrameAlgorithm(wrappers);
    declareSimpleAlgorithm(wrappers);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
