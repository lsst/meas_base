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
#include <pybind11/stl.h>

#include "lsst/pex/config/python.h"

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/PixelFlags.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(_pixelFlags) {
    py::module mod("_pixelFlags", "Python wrapper for afw _pixelFlags library");

    /* Module level */
    py::class_<PixelFlagsAlgorithm, std::shared_ptr<PixelFlagsAlgorithm>, SimpleAlgorithm> clsPixelFlagsAlgorithm(mod, "PixelFlagsAlgorithm");
    py::class_<PixelFlagsControl> clsPixelFlagsControl(mod, "PixelFlagsControl");

    /* Member types and enums */

    /* Constructors */
    clsPixelFlagsAlgorithm.def(py::init<PixelFlagsAlgorithm::Control const &,
                                        std::string const &,
                                        afw::table::Schema &>(),
                               "ctrl"_a, "name"_a, "schema"_a);

    clsPixelFlagsControl.def(py::init<>());

    /* Operators */

    /* Members */
    clsPixelFlagsAlgorithm.def("measure", &PixelFlagsAlgorithm::measure,
        "measRecord"_a, "exposure"_a);

    LSST_DECLARE_CONTROL_FIELD(clsPixelFlagsControl, PixelFlagsControl, masksFpAnywhere);
    LSST_DECLARE_CONTROL_FIELD(clsPixelFlagsControl, PixelFlagsControl, masksFpCenter);

    return mod.ptr();
}

}}}     // lsst::meas::base
