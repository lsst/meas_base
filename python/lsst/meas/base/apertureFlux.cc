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
#include "pybind11/stl.h"

#include <memory>

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/FluxUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyFluxAlgorithm =
        py::class_<ApertureFluxAlgorithm, std::shared_ptr<ApertureFluxAlgorithm>, SimpleAlgorithm>;
using PyFluxControl = py::class_<ApertureFluxControl>;
using PyFluxResult = py::class_<ApertureFluxResult, std::shared_ptr<ApertureFluxResult>, FluxResult>;
using PyFluxTransform =
        py::class_<ApertureFluxTransform, std::shared_ptr<ApertureFluxTransform>, BaseTransform>;

PyFluxControl declareFluxControl(py::module &mod) {
    PyFluxControl cls(mod, "ApertureFluxControl");

    LSST_DECLARE_CONTROL_FIELD(cls, ApertureFluxControl, radii);
    LSST_DECLARE_CONTROL_FIELD(cls, ApertureFluxControl, maxSincRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, ApertureFluxControl, shiftKernel);

    cls.def(py::init<>());

    return cls;
}

template <typename Image, class PyClass>
void declareComputeFluxes(PyClass &cls) {
    namespace afwEllipses = lsst::afw::geom::ellipses;
    using Control = ApertureFluxAlgorithm::Control;
    using Result = ApertureFluxAlgorithm::Result;

    cls.def_static("computeSincFlux",
                   (Result(*)(Image const &, afwEllipses::Ellipse const &, Control const &)) &
                           ApertureFluxAlgorithm::computeSincFlux,
                   "image"_a, "ellipse"_a, "ctrl"_a = Control());
    cls.def_static("computeNaiveFlux",
                   (Result(*)(Image const &, afwEllipses::Ellipse const &, Control const &)) &
                           ApertureFluxAlgorithm::computeNaiveFlux,
                   "image"_a, "ellipse"_a, "ctrl"_a = Control());
    cls.def_static("computeFlux", (Result(*)(Image const &, afwEllipses::Ellipse const &, Control const &)) &
                                          ApertureFluxAlgorithm::computeFlux,
                   "image"_a, "ellipse"_a, "ctrl"_a = Control());
}

PyFluxAlgorithm declareFluxAlgorithm(py::module &mod) {
    PyFluxAlgorithm cls(mod, "ApertureFluxAlgorithm");

    cls.attr("FAILURE") = py::cast(ApertureFluxAlgorithm::FAILURE);
    cls.attr("APERTURE_TRUNCATED") = py::cast(ApertureFluxAlgorithm::APERTURE_TRUNCATED);
    cls.attr("SINC_COEFFS_TRUNCATED") = py::cast(ApertureFluxAlgorithm::SINC_COEFFS_TRUNCATED);

    // constructor not wrapped because class is abstract

    declareComputeFluxes<lsst::afw::image::Image<double>>(cls);
    declareComputeFluxes<lsst::afw::image::MaskedImage<double>>(cls);
    declareComputeFluxes<lsst::afw::image::Image<float>>(cls);
    declareComputeFluxes<lsst::afw::image::MaskedImage<float>>(cls);

    cls.def("measure", &ApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &ApertureFluxAlgorithm::fail, "measRecord"_a, "error"_a = nullptr);
    cls.def_static("makeFieldPrefix", &ApertureFluxAlgorithm::makeFieldPrefix, "name"_a, "radius"_a);

    return cls;
}

void declareFluxResult(py::module &mod) {
    PyFluxResult cls(mod, "ApertureFluxResult");

    cls.def("getFlag", (bool (ApertureFluxResult::*)(unsigned int) const) & ApertureFluxResult::getFlag,
            "bit"_a);
    cls.def("getFlag",
            (bool (ApertureFluxResult::*)(std::string const &name) const) & ApertureFluxResult::getFlag,
            "name"_a);
    cls.def("setFlag", &ApertureFluxResult::setFlag, "index"_a, "value"_a);
    cls.def("unsetFlag", &ApertureFluxResult::unsetFlag, "index"_a);
}

PyFluxTransform declareFluxTransform(py::module &mod) {
    PyFluxTransform cls(mod, "ApertureFluxTransform");

    cls.def(py::init<ApertureFluxTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    cls.def("__call__", &ApertureFluxTransform::operator(), "inputCatalog"_a, "outputCatalog"_a, "wcs"_a,
            "calib"_a);

    return cls;
}

}  // <anonymous>

PYBIND11_PLUGIN(apertureFlux) {
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.fluxUtilities");
    py::module::import("lsst.meas.base.transform");

    py::module mod("apertureFlux");

    auto clsFluxControl = declareFluxControl(mod);
    auto clsFluxAlgorithm = declareFluxAlgorithm(mod);
    declareFluxResult(mod);
    auto clsFluxTransform = declareFluxTransform(mod);

    clsFluxAlgorithm.attr("Control") = clsFluxControl;
    // no need to make ApertureFluxControl::Result visible to Python
    clsFluxTransform.attr("Control") = clsFluxControl;

    python::declareAlgorithm<ApertureFluxAlgorithm, ApertureFluxControl, ApertureFluxTransform>(
            clsFluxAlgorithm, clsFluxControl, clsFluxTransform);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
