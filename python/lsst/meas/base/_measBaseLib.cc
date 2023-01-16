/*
 * LSST Data Management System
 * Copyright 2023 AURA/LSST.
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

namespace lsst {
namespace meas {
namespace base {
void wrapAlgorithm(lsst::utils::python::WrapperCollection &);
void wrapBlendedness(lsst::utils::python::WrapperCollection &);
void wrapCentroidUtilities(lsst::utils::python::WrapperCollection &);
void wrapCircularApertureFlux(lsst::utils::python::WrapperCollection &);
void wrapFlagHandler(lsst::utils::python::WrapperCollection &);
void wrapFluxUtilities(lsst::utils::python::WrapperCollection &);
void wrapTransform(lsst::utils::python::WrapperCollection &);
void wrapApertureFlux(lsst::utils::python::WrapperCollection &);
void wrapExceptions(lsst::utils::python::WrapperCollection &);
void wrapSincCoeffs(lsst::utils::python::WrapperCollection &);
void wrapShapeUtilities(lsst::utils::python::WrapperCollection &);
void wrapSdssShape(lsst::utils::python::WrapperCollection &);
void wrapSdssCentroid(lsst::utils::python::WrapperCollection &);
void wrapScaledApertureFlux(lsst::utils::python::WrapperCollection &);
void wrapPsfFlux(lsst::utils::python::WrapperCollection &);
void wrapNaiveCentroid(lsst::utils::python::WrapperCollection &);
void wrapLocalBackground(lsst::utils::python::WrapperCollection &);
void wrapGaussianFlux(lsst::utils::python::WrapperCollection &);
void wrapPeakLikelihoodFlux(lsst::utils::python::WrapperCollection &);
void wrapInputUtilities(lsst::utils::python::WrapperCollection &);
void wrapPixelFlags(lsst::utils::python::WrapperCollection &);

PYBIND11_MODULE(_measBaseLib, mod) {
    lsst::utils::python::WrapperCollection wrappers(mod, "lsst.meas.base");
    wrapExceptions(wrappers);  // No meas_base deps
    wrapSincCoeffs(wrappers);  // No meas_base deps
    wrapAlgorithm(wrappers);  // No meas_base deps
    wrapFlagHandler(wrappers);  // No meas_base deps
    wrapFluxUtilities(wrappers);  // No meas_base deps
    wrapTransform(wrappers);  // No meas_base deps
    wrapShapeUtilities(wrappers);   // No meas_base deps
    wrapInputUtilities(wrappers);  // No meas_base deps
    wrapApertureFlux(wrappers);  // Depends on algorithm, flagHandler, fluxUtilities, transform
    wrapScaledApertureFlux(wrappers);  // Depends on algorithm, fluxUtilities, transform
    wrapPsfFlux(wrappers);  // Depends on algorithm, fluxUtilities, transform, flagHandler
    wrapGaussianFlux(wrappers);  // Depends on algorithm, flagHandler, transform
    wrapPeakLikelihoodFlux(wrappers);  // Depends on algorithm, flagHandler, transform
    wrapBlendedness(wrappers);  // Depends on algorithm, flagHandler
    wrapPixelFlags(wrappers);  // Depends on algorithm
    wrapCentroidUtilities(wrappers);  // Depends on transform
    wrapCircularApertureFlux(wrappers);  // Depends on algorithm, apertureFlux
    wrapSdssCentroid(wrappers);  // Depends on algorithm, flagHandler, transform
    wrapNaiveCentroid(wrappers);  // Depends on algorithm, flagHandler, transform
    wrapSdssShape(wrappers);  // Depends on algorithm, flagHandler, centroidUtilities, fluxUtilities, shapeUtilities, transform
    wrapLocalBackground(wrappers);  // Depends on algorithm, flagHandler, transform, fluxUtilities
    wrappers.finish();
}
}  // namespace base
}  // namespace meas
}  // namespace lsst
