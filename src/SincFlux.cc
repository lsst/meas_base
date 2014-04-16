// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include "ndarray/eigen.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/meas/base/algorithms/SincFluxTemplates.h"
#include "lsst/meas/base/algorithms/detail/SincPhotometry.h"
#include "lsst/meas/base/SincFlux.h"

// Doxygen gets confused and generates warnings when trying to map the definitions here to their
// declarations, but we want to put the source code itself in the HTML docs, so we just tell it
// not to look for any documentation comments here.
/// @cond SOURCE_FILE

namespace lsst { namespace meas { namespace base {

SincFluxAlgorithm::ResultMapper SincFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    // calculate the needed coefficients 
    if (algorithms::photometry::fuzzyCompare<float>().isEqual(ctrl.ellipticity, 0.0)) {
        algorithms::photometry::SincCoeffs<float>::cache(ctrl.radius1, ctrl.radius2);
    }
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
SincFluxAlgorithm::Result SincFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Control const & ctrl
) {
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_PSF].doc,
            NO_PSF
        );
    }
    Result result;

    afw::geom::ellipses::Axes const axes(ctrl.radius2, ctrl.radius2*(1.0 - ctrl.ellipticity), ctrl.angle);
    std::pair<double, double> fluxes =
        algorithms::photometry::calculateSincApertureFlux(exposure.getMaskedImage(),
                                              afw::geom::ellipses::Ellipse(axes, center),
                                              ctrl.radius1/ctrl.radius2);
    double flux = fluxes.first;
    double fluxErr = fluxes.second;
    result.flux = flux;
    result.fluxSigma = fluxErr;

    //  End of meas_algorithms code
    return result;
}

template <typename T>
SincFluxAlgorithm::Result SincFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template SincFluxAlgorithm::Result SincFluxAlgorithm::apply(        \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    SincFluxAlgorithm::Result SincFluxAlgorithm::apply(                 \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    );                                                                  \
    template lsst::afw::image::Image<T>::Ptr algorithms::detail::calcImageRealSpace<T>(double const, double const,   \
                                                                           double const);                            \
    template lsst::afw::image::Image<T>::Ptr algorithms::detail::calcImageKSpaceReal<T>(double const, double const); \
    template lsst::afw::image::Image<T>::Ptr algorithms::detail::calcImageKSpaceCplx<T>(double const, double const,  \
                                                                            double const, double const);             \
    template std::pair<double, double>                                                                               \
    algorithms::photometry::calculateSincApertureFlux<lsst::afw::image::MaskedImage<T> >(                            \
        lsst::afw::image::MaskedImage<T> const&, lsst::afw::geom::ellipses::Ellipse const&, double const);           \
    template class algorithms::photometry::SincCoeffs<T>;

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

/// @endcond
