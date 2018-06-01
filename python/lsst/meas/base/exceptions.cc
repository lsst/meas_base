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

#include "lsst/meas/base/exceptions.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/exceptions/python/Exception.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(exceptions) {
    using pex::exceptions::python::declareException;
    using pex::exceptions::DomainError;
    using pex::exceptions::RuntimeError;

    py::module::import("lsst.pex.exceptions");

    py::module mod("exceptions");

    /* Module level */
    auto clsFatalAlgorithmError =
            declareException<FatalAlgorithmError, RuntimeError>(mod, "FatalAlgorithmError", "RuntimeError");
    auto clsMeasurementError =
            declareException<MeasurementError, RuntimeError>(mod, "MeasurementError", "RuntimeError");
    auto clsPixelValueError =
            declareException<PixelValueError, DomainError>(mod, "PixelValueError", "DomainError");

    /* Constructors */
    clsMeasurementError.def(py::init<std::string const &, std::size_t>(), "message"_a, "flagBit"_a);
    clsFatalAlgorithmError.def(py::init<std::string const &>(), "message"_a);
    clsPixelValueError.def(py::init<std::string const &>(), "message"_a);

    /* Members */
    clsMeasurementError.def("getFlagBit", &MeasurementError::getFlagBit);

    return mod.ptr();
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
