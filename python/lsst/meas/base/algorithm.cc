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

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/Algorithm.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_MODULE(algorithm, mod) {
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.table");

    py::class_<BaseAlgorithm, std::shared_ptr<BaseAlgorithm>> clsBaseAlgorithm(mod, "BaseAlgorithm");
    py::class_<SingleFrameAlgorithm, std::shared_ptr<SingleFrameAlgorithm>, BaseAlgorithm>
            clsSingleFrameAlgorithm(mod, "SingleFrameAlgorithm");
    py::class_<SimpleAlgorithm, std::shared_ptr<SimpleAlgorithm>, SingleFrameAlgorithm> clsSimpleAlgorithm(
            mod, "SimpleAlgorithm", py::multiple_inheritance());

    clsBaseAlgorithm.def("fail", &BaseAlgorithm::fail, "measRecord"_a, "error"_a = NULL);
    clsBaseAlgorithm.def("getLogName", &SimpleAlgorithm::getLogName);

    clsSingleFrameAlgorithm.def("measure", &SingleFrameAlgorithm::measure, "record"_a, "exposure"_a);

    clsSimpleAlgorithm.def("measureForced", &SimpleAlgorithm::measureForced, "measRecord"_a, "exposure"_a,
                           "refRecord"_a, "refWcs"_a);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
