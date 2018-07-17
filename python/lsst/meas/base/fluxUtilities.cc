/*
 * LSST Data Management System
 * Copyright 2008-2017  AURA/LSST.
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

#include "lsst/meas/base/ApertureFlux.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxResult = py::class_<FluxResult, std::shared_ptr<FluxResult>>;
using PyFluxResultKey = py::class_<FluxResultKey, std::shared_ptr<FluxResultKey>>;
using PyMagResult = py::class_<MagResult, std::shared_ptr<MagResult>>;
using PyMagResultKey = py::class_<MagResultKey, std::shared_ptr<MagResultKey>>;

void declareFluxResult(py::module &mod) {
    PyFluxResult cls(mod, "FluxResult");

    cls.def_readwrite("flux", &FluxResult::flux);
    cls.def_readwrite("fluxSigma", &FluxResult::fluxSigma);
}

void declareFluxResultKey(py::module &mod) {
    PyFluxResultKey cls(mod, "FluxResultKey");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::Key<meas::base::Flux> const &, afw::table::Key<FluxErrElement> const &>(),
            "flux"_a, "fluxSigma"_a);
    cls.def(py::init<afw::table::SubSchema const &>());

    cls.def("__eq__", &FluxResultKey::operator==, py::is_operator());
    cls.def("__ne__", &FluxResultKey::operator!=, py::is_operator());

    cls.def("get", &FluxResultKey::get);
    cls.def("set", &FluxResultKey::set);
    cls.def_static("addFields", &FluxResultKey::addFields, "schema"_a, "name"_a, "doc"_a);
    cls.def("isValid", &FluxResultKey::isValid);
    cls.def("getFlux", &FluxResultKey::getFlux);
    cls.def("getFluxSigma", &FluxResultKey::getFluxSigma);
}

void declareMagResult(py::module &mod) {
    PyMagResult cls(mod, "MagResult");

    cls.def_readwrite("mag", &MagResult::mag);
    cls.def_readwrite("magErr", &MagResult::magErr);
}

void declareMagResultKey(py::module &mod) {
    PyMagResultKey cls(mod, "MagResultKey");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::SubSchema const &>());

    cls.def("get", &MagResultKey::get);
    cls.def("set",
            (void (MagResultKey::*)(afw::table::BaseRecord &, MagResult const &) const) & MagResultKey::set);
    cls.def("set",
            (void (MagResultKey::*)(afw::table::BaseRecord &, std::pair<double, double> const &) const) &
                    MagResultKey::set);
    cls.def_static("addFields", &MagResultKey::addFields, "schema"_a, "name"_a);
}

}  // namespace

PYBIND11_PLUGIN(fluxUtilities) {
    py::module::import("lsst.afw.table");

    py::module mod("fluxUtilities");

    declareFluxResult(mod);
    declareFluxResultKey(mod);
    declareMagResult(mod);
    declareMagResultKey(mod);

    return mod.ptr();
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
