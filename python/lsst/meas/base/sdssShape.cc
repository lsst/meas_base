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

#include "lsst/pex/config/pybind11.h"
#include "lsst/meas/base/pybind11.h"

#include "lsst/meas/base/SdssShape.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(_sdssShape) {
    py::module mod("_sdssShape", "Python wrapper for afw _sdssShape library");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    /* Module level */
    py::class_<SdssShapeResult, CentroidResult, FluxResult, ShapeResult> clsSdssShapeResult(mod, "SdssShapeResult");
    py::class_<SdssShapeResultKey> clsSdssShapeResultKey(mod, "SdssShapeResultKey");
    py::class_<SdssShapeAlgorithm, SimpleAlgorithm> clsSdssShapeAlgorithm(mod, "SdssShapeAlgorithm");
    py::class_<SdssShapeControl> clsSdssShapeControl(mod, "SdssShapeControl");
    py::class_<SdssShapeTransform> clsSdssShapeTransform(mod, "SdssShapeTransform");

    mod.def("makeShapeTransformMatrix", &makeShapeTransformMatrix, "xform"_a);

    /* Member types and enums */
    // Anonymous enum cannot be represented directly
    clsSdssShapeAlgorithm.attr("FAILURE") = py::cast(static_cast<int>(SdssShapeAlgorithm::FAILURE));
    clsSdssShapeAlgorithm.attr("UNWEIGHTED_BAD") = py::cast(static_cast<int>(SdssShapeAlgorithm::UNWEIGHTED_BAD));
    clsSdssShapeAlgorithm.attr("UNWEIGHTED") = py::cast(static_cast<int>(SdssShapeAlgorithm::UNWEIGHTED));
    clsSdssShapeAlgorithm.attr("SHIFT") = py::cast(static_cast<int>(SdssShapeAlgorithm::SHIFT));
    clsSdssShapeAlgorithm.attr("MAXITER") = py::cast(static_cast<int>(SdssShapeAlgorithm::MAXITER));
    clsSdssShapeAlgorithm.attr("PSF_SHAPE_BAD") = py::cast(static_cast<int>(SdssShapeAlgorithm::PSF_SHAPE_BAD));
    clsSdssShapeAlgorithm.attr("N_FLAGS") = py::cast(static_cast<int>(SdssShapeAlgorithm::N_FLAGS));

    /* Constructors */
    clsSdssShapeResultKey.def(py::init<afw::table::SubSchema const &>(),
                              "s"_a);

    /* Operators */

    /* Members */
    python::declareAlgorithm<SdssShapeAlgorithm,
                             SdssShapeControl,
                             SdssShapeTransform>(clsSdssShapeAlgorithm,
                                                 clsSdssShapeControl,
                                                 clsSdssShapeTransform);

    LSST_DECLARE_CONTROL_FIELD(clsSdssShapeControl, SdssShapeControl, background);
    LSST_DECLARE_CONTROL_FIELD(clsSdssShapeControl, SdssShapeControl, doMeasurePsf);
    LSST_DECLARE_CONTROL_FIELD(clsSdssShapeControl, SdssShapeControl, maxIter);
    LSST_DECLARE_CONTROL_FIELD(clsSdssShapeControl, SdssShapeControl, maxShift);
    LSST_DECLARE_CONTROL_FIELD(clsSdssShapeControl, SdssShapeControl, tol1);
    LSST_DECLARE_CONTROL_FIELD(clsSdssShapeControl, SdssShapeControl, tol2);

    // TODO: pybind11 documentation for this method says it's a workaround for Swig
    clsSdssShapeResult.def("getFlag", &SdssShapeResult::getFlag);

    clsSdssShapeResult.def_readwrite("flux_xx_Cov", &SdssShapeResult::flux_xx_Cov);
    clsSdssShapeResult.def_readwrite("flux_yy_Cov", &SdssShapeResult::flux_yy_Cov);
    clsSdssShapeResult.def_readwrite("flux_xy_Cov", &SdssShapeResult::flux_xy_Cov);

    clsSdssShapeResultKey.def("get", &SdssShapeResultKey::get, "record"_a);
    clsSdssShapeResultKey.def("getPsfShape",
                              (afwGeom::ellipses::Quadrupole (SdssShapeResultKey::*)
                                    (afw::table::BaseRecord const &)) &SdssShapeResultKey::getPsfShape,
                              "record"_a);

    return mod.ptr();
}

}}}     // lsst::meas::base
