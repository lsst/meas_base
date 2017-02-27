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

#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"
#include "ndarray/converter.h"

#include "lsst/meas/base/CentroidUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(centroidUtilities) {
    py::module mod("centroidUtilities", "Python wrapper for afw _centroidUtilities library");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    /* Module level */
    py::class_<CentroidChecker> clsCentroidChecker(mod, "CentroidChecker");
    py::class_<CentroidResult> clsCentroidResult(mod, "CentroidResult");
    py::class_<CentroidResultKey> clsCentroidResultKey(mod, "CentroidResultKey");

    /* Member types and enums */

    /* Constructors */
    clsCentroidChecker.def(py::init<afw::table::Schema &,
                                    std::string const &,
                                    bool,
                                    double>(),
                           "schema"_a, "name"_a, "inside"_a=true, "maxDistFromPeak"_a=-1.0);

    clsCentroidResultKey.def(py::init<afw::table::SubSchema const &>(), "s"_a);

    /* Operators */
    clsCentroidChecker.def("__call__", [](CentroidChecker const & self,
                                          afw::table::SourceRecord & record) {
            return self(record);
        }, "record"_a);

    /* Members */
    clsCentroidResult.def("getCentroid", &CentroidResult::getCentroid);
    clsCentroidResult.def("getCentroidErr", &CentroidResult::getCentroidErr);

    clsCentroidResult.def_readwrite("x", &CentroidResult::x);
    clsCentroidResult.def_readwrite("y", &CentroidResult::y);

    clsCentroidResultKey.def("get", &CentroidResultKey::get, "record"_a);
    clsCentroidResultKey.def("getCentroidErr", &CentroidResultKey::getCentroidErr);

    return mod.ptr();
}

}}}     // lsst::meas::base
