/*
 * This file is part of meas_base.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

using cpputils::python::WrapperCollection;

void wrapFluxUtilities(WrapperCollection&);
void wrapAlgorithm(WrapperCollection&);
void wrapApertureFlow(WrapperCollection&);
void wrapBlendedness(WrapperCollection&);
void wrapCentroidUtilities(WrapperCollection&);
void wrapCircularApertureFlux(WrapperCollection&);
void wrapExceptions(WrapperCollection&);
void wrapFlagHandler(WrapperCollection&);
void wrapGaussianFlux(WrapperCollection &);
void wrapInputUtilities(WrapperCollection&);
void wrapLocalBackground(WrapperCollection&);
void wrapPeakLikelihoodFlux(WrapperCollection&);
void wrapPixelFLags(WrapperCollection&);
void wrapPsfFlux(WrapperCollection&);
void wrapScaledApertureFlux(WrapperCollection&);
void wrapSddsCentroid(WrapperCollection&);
void wrapShapeUtilities(WrapperCollection&);
void wrapSincCoeffs(WrapperCollection&);
void wrapSsdsShape(WrapperCollection&);
void wrapTransform(WrapperCollection&);
void wrapCalcCompensatedGaussian(WrapperCollection&);

PYBIND11_MODULE(_measBaseLib, mod) {
    lsst::cpputils::python::WrapperCollection wrappers(mod, "lsst.meas.base");

    wrappers.addInheritanceDependency("lsst.afw.geom");
    wrappers.addInheritanceDependency("lsst.afw.image");
    wrappers.addInheritanceDependency("lsst.afw.table");
    wrappers.addInheritanceDependency("lsst.pex.exceptions");

    wrappers.addSignatureDependency("lsst.daf.base");

    wrapExceptions(wrappers);
    wrapFlagHandler(wrappers);
    wrapFluxUtilities(wrappers);
    wrapTransform(wrappers);
    wrapAlgorithm(wrappers);
    wrapApertureFlow(wrappers);
    wrapBlendedness(wrappers);
    wrapCentroidUtilities(wrappers);
    wrapCircularApertureFlux(wrappers);
    wrapGaussianFlux(wrappers);
    wrapInputUtilities(wrappers);
    wrapLocalBackground(wrappers);
    wrapPeakLikelihoodFlux(wrappers);
    wrapPixelFLags(wrappers);
    wrapPsfFlux(wrappers);
    wrapScaledApertureFlux(wrappers);
    wrapSddsCentroid(wrappers);
    wrapShapeUtilities(wrappers);
    wrapSincCoeffs(wrappers);
    wrapSsdsShape(wrappers);
    wrapCalcCompensatedGaussian(wrappers);
    wrappers.finish();
}
}  // namespace base
}  // namespace meas
}  // namespace lsst
