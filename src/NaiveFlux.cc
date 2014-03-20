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
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/base/algorithms/NaiveFluxTemplates.h"
#include "lsst/meas/base/NaiveFlux.h"

// Doxygen gets confused and generates warnings when trying to map the definitions here to their
// declarations, but we want to put the source code itself in the HTML docs, so we just tell it
// not to look for any documentation comments here.
/// @cond SOURCE_FILE

namespace lsst { namespace meas { namespace base {

NaiveFluxAlgorithm::ResultMapper NaiveFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
NaiveFluxAlgorithm::Result NaiveFluxAlgorithm::apply(
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
    typename afw::image::Exposure<T>::MaskedImageT const& mimage = exposure.getMaskedImage();

    double const xcen = center.getX();   ///< object's column position
    double const ycen = center.getY();   ///< object's row position

    int const ixcen = afw::image::positionToIndex(xcen);
    int const iycen = afw::image::positionToIndex(ycen);

    // BBox for data image
    afw::geom::BoxI imageBBox(mimage.getBBox(afw::image::PARENT));

    /* ******************************************************* */
    // Aperture flux
    algorithms::FootprintFlux<typename afw::image::Exposure<T>::MaskedImageT> fluxFunctor(mimage);
    afw::detection::Footprint const foot(
        afw::geom::PointI(ixcen, iycen), 
        ctrl.radius, 
        imageBBox
        );
    fluxFunctor.apply(foot);

    result.flux = fluxFunctor.getSum();
    result.fluxSigma = ::sqrt(fluxFunctor.getSumVar());
    return result;
}

template <typename T>
NaiveFluxAlgorithm::Result NaiveFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template NaiveFluxAlgorithm::Result NaiveFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    NaiveFluxAlgorithm::Result NaiveFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

/// @endcond
