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

#include "lsst/meas/base/FlagHandler.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(flagHandler) {
    py::module mod("flagHandler", "Python wrapper for afw _flagHandler library");

    /* Module level */
    py::class_<FlagDefinition> clsFlagDefinition(mod, "FlagDefinition");
    py::class_<FlagHandler> clsFlagHandler(mod, "FlagHandler");

    /* Member types and enums */
    // Anonymous enum cannot be represented directly
    clsFlagHandler.attr("FAILURE") = py::cast(static_cast<int>(FlagHandler::FAILURE));

    /* Constructors */
    clsFlagDefinition.def(py::init<>());
    clsFlagDefinition.def(py::init<std::string, std::string>(), "_name"_a, "_doc"_a);

    /* Operators */

    /* Members */
    clsFlagDefinition.def_readwrite("doc", &FlagDefinition::doc);
    clsFlagDefinition.def_readwrite("name", &FlagDefinition::name);

    clsFlagHandler.def("getDefinition", &FlagHandler::getDefinition, "i"_a);
    clsFlagHandler.def("getValue", &FlagHandler::getValue, "record"_a, "i"_a);
    // Default value actually NULL, not nullptr, but python doesn't like passing integer as an exception
    clsFlagHandler.def("handleFailure", &FlagHandler::handleFailure, "record"_a, "error"_a=nullptr);
    clsFlagHandler.def_static("addFields",
                              (FlagHandler (*)(afw::table::Schema &,
                                               std::string const &,
                                               std::vector<FlagDefinition> const *)) &FlagHandler::addFields,
                              "schema"_a, "prefix"_a, "flagDefs"_a);

    return mod.ptr();
}

}}}     // lsst::meas::base
