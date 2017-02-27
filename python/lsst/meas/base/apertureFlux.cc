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

#include "lsst/meas/base/ApertureFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {
    template <typename Image, class PyClass>
    void declareComputeFluxes(PyClass & cls) {
        namespace afwEllipses = lsst::afw::geom::ellipses;
        using Control = ApertureFluxAlgorithm::Control;
        using Result = ApertureFluxAlgorithm::Result;

        cls.def_static("computeSincFlux", (Result (*) (Image const &,
                                                       afwEllipses::Ellipse const &,
                                                       Control const &)) &ApertureFluxAlgorithm::computeSincFlux,
            "image"_a, "ellipse"_a, "ctrl"_a=Control());
        cls.def_static("computeNaiveFlux", (Result (*) (Image const &,
                                                        afwEllipses::Ellipse const &,
                                                        Control const &)) &ApertureFluxAlgorithm::computeNaiveFlux,
            "image"_a, "ellipse"_a, "ctrl"_a=Control());
        cls.def_static("computeFlux", (Result (*) (Image const &,
                                                   afwEllipses::Ellipse const &,
                                                   Control const &)) &ApertureFluxAlgorithm::computeFlux,
            "image"_a, "ellipse"_a, "ctrl"_a=Control());
    }
}

PYBIND11_PLUGIN(apertureFlux) {
    py::module mod("apertureFlux");

    /* Module level */
    py::class_<ApertureFluxAlgorithm, std::shared_ptr<ApertureFluxAlgorithm>, SimpleAlgorithm> clsApertureFluxAlgorithm(mod, "ApertureFluxAlgorithm");
    py::class_<ApertureFluxControl> clsApertureFluxControl(mod, "ApertureFluxControl");
    py::class_<ApertureFluxResult, FluxResult> clsApertureFluxResult(mod, "ApertureFluxResult");
    py::class_<ApertureFluxTransform> clsApertureFluxTransform(mod, "ApertureFluxTransform");

    /* Member types and enums */
    clsApertureFluxAlgorithm.attr("Control") = clsApertureFluxControl;

    // ApertureFluxAlgorithm::Control wrapped in apertureFlux.py
    py::enum_<ApertureFluxAlgorithm::FlagBits>(clsApertureFluxAlgorithm, "FlagBits")
        .value("FAILURE", ApertureFluxAlgorithm::FlagBits::FAILURE)
        .value("APERTURE_TRUNCATED", ApertureFluxAlgorithm::FlagBits::APERTURE_TRUNCATED)
        .value("SINC_COEFFS_TRUNCATED", ApertureFluxAlgorithm::FlagBits::SINC_COEFFS_TRUNCATED)
        .value("N_FLAGS", ApertureFluxAlgorithm::FlagBits::N_FLAGS)
        .export_values();

    /* Members */
    python::declareAlgorithm<ApertureFluxAlgorithm,
                             ApertureFluxControl,
                             ApertureFluxTransform>(clsApertureFluxAlgorithm,
                                                    clsApertureFluxControl,
                                                    clsApertureFluxTransform);

    declareComputeFluxes<lsst::afw::image::Image<double>>(clsApertureFluxAlgorithm);
    declareComputeFluxes<lsst::afw::image::MaskedImage<double>>(clsApertureFluxAlgorithm);
    declareComputeFluxes<lsst::afw::image::Image<float>>(clsApertureFluxAlgorithm);
    declareComputeFluxes<lsst::afw::image::MaskedImage<float>>(clsApertureFluxAlgorithm);

    clsApertureFluxAlgorithm.def_static("makeFieldPrefix", &ApertureFluxAlgorithm::makeFieldPrefix, 
        "name"_a, "radius"_a);
    
    LSST_DECLARE_CONTROL_FIELD(clsApertureFluxControl, ApertureFluxControl, maxSincRadius);
    LSST_DECLARE_CONTROL_FIELD(clsApertureFluxControl, ApertureFluxControl, radii);
    LSST_DECLARE_CONTROL_FIELD(clsApertureFluxControl, ApertureFluxControl, shiftKernel);
    
    clsApertureFluxResult.def("getFlag", &ApertureFluxResult::getFlag, "bit"_a);

    return mod.ptr();
}

}}}     // lsst::meas::base
