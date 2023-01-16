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
#include "lsst/utils/python.h"

#include "lsst/meas/base/SincCoeffs.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

template <typename T>
void declareSincCoeffs(lsst::utils::python::WrapperCollection &wrappers, std::string const& suffix) {
    std::string name = "SincCoeffs" + suffix;

    wrappers.wrapType(
            py::class_<SincCoeffs<T>>(wrappers.module, name.c_str()),
            [](auto &mod, auto &cls) {
                cls.def_static("cache", &SincCoeffs<T>::cache, "rInner"_a, "rOuter"_a);
                cls.def_static("get", &SincCoeffs<T>::get, "outerEllipse"_a, "innerRadiusFactor"_a);
    });
}

}  // namespace

void wrapSincCoeffs(lsst::utils::python::WrapperCollection &wrappers) {
    wrappers.addInheritanceDependency("lsst.afw.geom");
    wrappers.addInheritanceDependency("lsst.afw.image");

    declareSincCoeffs<float>(wrappers, "F");
    declareSincCoeffs<double>(wrappers, "D");
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
