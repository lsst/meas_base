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

#include "pybind11/pybind11.h"

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/ShapeUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(shapeUtilities) {
    py::module mod("shapeUtilities");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    /* Module level */
    py::class_<ShapeResult> clsShapeResult(mod, "ShapeResult");
    py::class_<ShapeResultKey> clsShapeResultKey(mod, "ShapeResultKey");

    mod.def("makeShapeTransformMatrix", &makeShapeTransformMatrix, "xform"_a);

    /* Constructors */
    clsShapeResultKey.def(py::init<afw::table::SubSchema const &>(), "s"_a);

    /* Members */
    clsShapeResult.def("getShape", &ShapeResult::getShape);
    clsShapeResult.def("getShapeErr",
        (lsst::meas::base::ShapeCov const (ShapeResult::*)() const) &ShapeResult::getShapeErr);

    clsShapeResult.def_readwrite("xx", &ShapeResult::xx);
    clsShapeResult.def_readwrite("yy", &ShapeResult::yy);
    clsShapeResult.def_readwrite("xy", &ShapeResult::xy);
    clsShapeResult.def_readwrite("xxSigma", &ShapeResult::xxSigma);
    clsShapeResult.def_readwrite("yySigma", &ShapeResult::yySigma);
    clsShapeResult.def_readwrite("xySigma", &ShapeResult::xySigma);
    clsShapeResult.def_readwrite("xx_yy_Cov", &ShapeResult::xx_yy_Cov);
    clsShapeResult.def_readwrite("xx_xy_Cov", &ShapeResult::xx_xy_Cov);
    clsShapeResult.def_readwrite("yy_xy_Cov", &ShapeResult::yy_xy_Cov);

    clsShapeResultKey.def("get", &ShapeResultKey::get, "record"_a);

    return mod.ptr();
}

}}}     // lsst::meas::base
