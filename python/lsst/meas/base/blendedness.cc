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

#include "lsst/meas/base/Blendedness.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(_blendedness) {
    py::module mod("_blendedness", "Python wrapper for afw _blendedness library");

    /* Module level */

    /* Member types and enums */
    py::class_<BlendednessAlgorithm, std::shared_ptr<BlendednessAlgorithm>, SimpleAlgorithm> clsBlendednessAlgorithm(mod, "BlendednessAlgorithm");
    py::class_<BlendednessControl> clsBlendednessControl(mod, "BlendednessControl");

    /* Constructors */
    clsBlendednessAlgorithm.def(py::init<BlendednessAlgorithm::Control const &,
                                         std::string const &,
                                         afw::table::Schema &>(),
                                "ctrl"_a, "name"_a, "schema"_a);

    clsBlendednessControl.def(py::init<>());

    /* Operators */

    /* Members */
    python::declareAlgorithm<BlendednessAlgorithm, BlendednessControl>(
        clsBlendednessAlgorithm, clsBlendednessControl);
    clsBlendednessAlgorithm.def("measureChildPixels", &BlendednessAlgorithm::measureChildPixels,
        "image"_a, "child"_a);
    clsBlendednessAlgorithm.def("measureParentPixels", &BlendednessAlgorithm::measureParentPixels,
        "image"_a, "child"_a);

    LSST_DECLARE_CONTROL_FIELD(clsBlendednessControl, BlendednessControl, doFlux);
    LSST_DECLARE_CONTROL_FIELD(clsBlendednessControl, BlendednessControl, doOld);
    LSST_DECLARE_CONTROL_FIELD(clsBlendednessControl, BlendednessControl, doShape);
    LSST_DECLARE_CONTROL_FIELD(clsBlendednessControl, BlendednessControl, nSigmaWeightMax);

    return mod.ptr();
}

}}}     // lsst::meas::base
