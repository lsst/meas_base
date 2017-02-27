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
#include "pybind11/stl.h"

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/PsfFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(psfFlux) {
    py::module mod("psfFlux");

    /* Module level */
    py::class_<PsfFluxAlgorithm, std::shared_ptr<PsfFluxAlgorithm>, SimpleAlgorithm> clsPsfFluxAlgorithm(mod, "PsfFluxAlgorithm");
    py::class_<PsfFluxControl> clsPsfFluxControl(mod, "PsfFluxControl");
    py::class_<PsfFluxTransform> clsPsfFluxTransform(mod, "PsfFluxTransform");

    /* Member types and enums */
    // Anonymous enum cannot be represented directly
    clsPsfFluxAlgorithm.attr("FAILURE") = py::cast(static_cast<int>(PsfFluxAlgorithm::FAILURE));
    clsPsfFluxAlgorithm.attr("NO_GOOD_PIXELS") = py::cast(static_cast<int>(PsfFluxAlgorithm::NO_GOOD_PIXELS));
    clsPsfFluxAlgorithm.attr("EDGE") = py::cast(static_cast<int>(PsfFluxAlgorithm::EDGE));
    clsPsfFluxAlgorithm.attr("N_FLAGS") = py::cast(static_cast<int>(PsfFluxAlgorithm::N_FLAGS));

    /* Members */
    python::declareAlgorithm<PsfFluxAlgorithm,
                             PsfFluxControl,
                             PsfFluxTransform>(clsPsfFluxAlgorithm,
                                               clsPsfFluxControl,
                                               clsPsfFluxTransform);

    LSST_DECLARE_CONTROL_FIELD(clsPsfFluxControl, PsfFluxControl, badMaskPlanes);

    return mod.ptr();
}

}}}     // lsst::meas::base
