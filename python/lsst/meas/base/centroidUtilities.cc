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

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/meas/base/CentroidUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyCentroidChecker = py::class_<CentroidChecker>;
using PyCentroidResult = py::class_<CentroidResult, std::shared_ptr<CentroidResult>>;
using PyCentroidResultKey = py::class_<CentroidResultKey>;
using PyCentroidTransform = py::class_<CentroidTransform, std::shared_ptr<CentroidTransform>, BaseTransform>;

void declareCentroidResult(py::module & mod) {
    PyCentroidResult cls(mod, "CentroidResult");

    cls.def_readwrite("x", &CentroidResult::x);
    cls.def_readwrite("y", &CentroidResult::y);
    cls.def_readwrite("xSigma", &CentroidResult::xSigma);
    cls.def_readwrite("ySigma", &CentroidResult::ySigma);
    cls.def_readwrite("x_y_Cov", &CentroidResult::x_y_Cov);

    cls.def(py::init<>());
    cls.def(py::init<CentroidElement, CentroidElement, CentroidCov const &>(), "x"_a, "y"_a, "matrix"_a);
    cls.def(py::init<CentroidElement, CentroidElement, ErrElement, ErrElement>(), "x"_a, "y"_a, "xSigma"_a,
            "ySigma"_a);

    cls.def("getCentroid", &CentroidResult::getCentroid);
    cls.def("setCentroid", &CentroidResult::setCentroid, "centroid"_a);
    cls.def("getPoint", &CentroidResult::getPoint);
    cls.def("getCentroidErr", &CentroidResult::getCentroidErr);
    cls.def("setCentroidErr",
            (void (CentroidResult::*)(CentroidCov const &)) & CentroidResult::setCentroidErr, "matrix"_a);
    cls.def("setCentroidErr",
            (void (CentroidResult::*)(ErrElement, ErrElement)) & CentroidResult::setCentroidErr, "xSigma"_a,
            "ySigma"_a);
}

void declareCentroidResultKey(py::module & mod) {
    PyCentroidResultKey cls(mod, "CentroidResultKey");

    cls.def(py::init<>());
    // TODO make this work or document it not being needed
    // cls.def(py::init<afw::table::PointKey<CentroidElement> const &,
    //                  afw::table::CovarianceMatrixKey<ErrElement, 2> const &>(),
    //         "centroid"_a);
    cls.def(py::init<afw::table::SubSchema const &>(), "subSchema"_a);

    cls.def("__eq__", &CentroidResultKey::operator==, py::is_operator());
    cls.def("__nq__", &CentroidResultKey::operator!=, py::is_operator());

    cls.def("get", &CentroidResultKey::get, "record"_a);
    cls.def("set", &CentroidResultKey::set, "record"_a, "value"_a);
    cls.def("isValid", &CentroidResultKey::isValid);
    cls.def("getCentroid", &CentroidResultKey::getCentroid);
    cls.def("getCentroidErr", &CentroidResultKey::getCentroidErr);
    cls.def("getX", &CentroidResultKey::getX);
    cls.def("getY", &CentroidResultKey::getY);
}

void declareCentroidTransform(py::module & mod) {
    PyCentroidTransform cls(mod, "CentroidTransform");

    cls.def(py::init<std::string const &, afw::table::SchemaMapper &>(), "name"_a, "mapper"_a);

    cls.def("__call__", &CentroidTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
            "calib"_a);
}

void declareCentroidChecker(py::module &mod) {
    PyCentroidChecker cls(mod, "CentroidChecker");

    cls.def(py::init<afw::table::Schema &, std::string const &, bool, double>(), "schema"_a, "name"_a,
            "inside"_a = true, "maxDistFromPeak"_a = -1.0);

    cls.def("__call__", &CentroidChecker::operator(), "record"_a);
}

}  // <anonymous>

PYBIND11_PLUGIN(centroidUtilities) {
    py::module::import("lsst.meas.base.transform");

    py::module mod("centroidUtilities");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declareCentroidResult(mod);
    declareCentroidResultKey(mod);
    declareCentroidTransform(mod);
    declareCentroidChecker(mod);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
