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

#include "lsst/meas/base/exceptions.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/exceptions/python/Exception.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

void wrapExceptions(lsst::utils::python::WrapperCollection &wrappers) {
    using pex::exceptions::DomainError;
    using pex::exceptions::RuntimeError;

    wrappers.addSignatureDependency("lsst.pex.exceptions");

    auto clsFatalAlgorithmError = wrappers.wrapException<FatalAlgorithmError, RuntimeError>("FatalAlgorithmError", "RuntimeError");
    auto clsMeasurementError = wrappers.wrapException<MeasurementError, RuntimeError>("MeasurementError", "RuntimeError");
    auto clsPixelValueError = wrappers.wrapException<PixelValueError, DomainError>("PixelValueError", "DomainError");

    clsMeasurementError.def(py::init<std::string const &, std::size_t>(), "message"_a, "flagBit"_a);
    clsFatalAlgorithmError.def(py::init<std::string const &>(), "message"_a);
    clsPixelValueError.def(py::init<std::string const &>(), "message"_a);

    clsMeasurementError.def("getFlagBit", &MeasurementError::getFlagBit);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
