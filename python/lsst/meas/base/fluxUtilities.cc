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

#include "lsst/meas/base/ApertureFlux.h"

#include "lsst/afw/table/python/functorKey.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace meas {
namespace base {

PYBIND11_PLUGIN(fluxUtilities) {
    py::module mod("fluxUtilities", "Python wrapper for afw _fluxUtilities library");

    py::class_<FluxResult> clsFluxResult(mod, "FluxResult");

    clsFluxResult.def_readwrite("flux", &FluxResult::flux);
    clsFluxResult.def_readwrite("fluxSigma", &FluxResult::fluxSigma);

    afw::table::python::declareFunctorKeys<FluxResult>(mod, "FluxResult");

    py::class_<FluxResultKey, std::shared_ptr<FluxResultKey>, afw::table::FunctorKey<FluxResult>> clsFluxResultKey(mod, "FluxResultKey");

    clsFluxResultKey.def(py::init<>());
    clsFluxResultKey.def(py::init<afw::table::Key<meas::base::Flux> const &, afw::table::Key<FluxErrElement> const &>(),
            "flux"_a, "fluxSigma"_a);
    clsFluxResultKey.def(py::init<afw::table::SubSchema const &>());

    clsFluxResultKey.def("__eq__", &FluxResultKey::operator==, py::is_operator());
    clsFluxResultKey.def("__ne__", &FluxResultKey::operator!=, py::is_operator());

    clsFluxResultKey.def("get", &FluxResultKey::get);
    clsFluxResultKey.def("set", &FluxResultKey::set);
    clsFluxResultKey.def_static("addFields", &FluxResultKey::addFields,
            "schema"_a, "name"_a, "doc"_a);
    clsFluxResultKey.def("isValid", &FluxResultKey::isValid);
    clsFluxResultKey.def("getFlux", &FluxResultKey::getFlux);
    clsFluxResultKey.def("getFluxSigma", &FluxResultKey::getFluxSigma);

    py::class_<MagResult> clsMagResult(mod, "MagResult");

    clsMagResult.def_readwrite("mag", &MagResult::mag);
    clsMagResult.def_readwrite("magErr", &MagResult::magErr);

    afw::table::python::declareFunctorKeys<MagResult>(mod, "MagResult");

    py::class_<MagResultKey, std::shared_ptr<MagResultKey>, afw::table::FunctorKey<MagResult>> clsMagResultKey(mod, "MagResultKey");

    clsMagResultKey.def(py::init<>());
    clsMagResultKey.def(py::init<afw::table::SubSchema const &>());

    clsMagResultKey.def("get", &MagResultKey::get);
    clsMagResultKey.def("set", (void (MagResultKey::*)(afw::table::BaseRecord &, MagResult const &) const) &MagResultKey::set);
    clsMagResultKey.def("set", (void (MagResultKey::*)(afw::table::BaseRecord &, std::pair<double,double> const &) const) &MagResultKey::set);
    clsMagResultKey.def_static("addFields", &MagResultKey::addFields,
            "schema"_a, "name"_a);

    return mod.ptr();
}

}}}     // lsst::meas::base
