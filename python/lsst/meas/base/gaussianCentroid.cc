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

#include "lsst/pex/config/pybind11.h"
#include "lsst/meas/base/pybind11.h"

#include "lsst/meas/base/GaussianCentroid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(_gaussianCentroid) {
    py::module mod("_gaussianCentroid", "Python wrapper for afw _gaussianCentroid library");

    /* Module level */
    py::class_<GaussianCentroidAlgorithm, SimpleAlgorithm> clsGaussianCentroidAlgorithm(mod, "GaussianCentroidAlgorithm");
    py::class_<GaussianCentroidControl> clsGaussianCentroidControl(mod, "GaussianCentroidControl");
    py::class_<GaussianCentroidTransform> clsGaussianCentroidTransform(mod, "GaussianCentroidTransform");

    /* Member types and enums */

    /* Constructors */

    /* Operators */


    /* Members */
    python::declareAlgorithm<GaussianCentroidAlgorithm,
                             GaussianCentroidControl,
                             GaussianCentroidTransform>(clsGaussianCentroidAlgorithm,
                                                        clsGaussianCentroidControl,
                                                        clsGaussianCentroidTransform);
    LSST_DECLARE_CONTROL_FIELD(clsGaussianCentroidControl, GaussianCentroidControl, doFootprintCheck);
    LSST_DECLARE_CONTROL_FIELD(clsGaussianCentroidControl, GaussianCentroidControl, maxDistToPeak);

    return mod.ptr();
}

}}}     // lsst::meas::base
