// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#include <numeric>

#include "ndarray/eigen.h"

#include "lsst/afw/math/offsetImage.h"
#include "lsst/afw/geom/ellipses/PixelRegion.h"
#include "lsst/meas/base/SincCoeffs.h"
#include "lsst/meas/base/ApertureFlux.h"

namespace lsst { namespace meas { namespace base {

ApertureFluxComponent::ApertureFluxComponent(int size) :
    flux(ndarray::allocate(size)),
    fluxSigma(ndarray::allocate(size))
{
    flux.deep() = std::numeric_limits<Flux>::quiet_NaN();
    fluxSigma.deep() = std::numeric_limits<FluxErrElement>::quiet_NaN();
}

ApertureFluxComponentMapper::ApertureFluxComponentMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    std::vector<double> const & radii
) :
    _flux(
        afw::table::ArrayKey<Flux>::addFields(
            schema,
            prefix + "_flux",
            "flux within %f pixel aperture",
            "dn",
            radii
        )
    ),
    _fluxSigma(
        afw::table::ArrayKey<Flux>::addFields(
            schema,
            prefix + "_fluxSigma",
            "1-sigma uncertainty on flux within %f pixel aperture",
            "dn",
            radii
        )
    )
{}

void ApertureFluxComponentMapper::apply(
    afw::table::BaseRecord & record,
    ApertureFluxComponent const & result
) const {
    record.set(_flux, result.flux);
    record.set(_fluxSigma, result.fluxSigma);
}

template <typename T>
Flux ApertureFluxAlgorithm::computeSincFlux(
    afw::image::Image<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    CONST_PTR(afw::image::Image<T>) cImage = SincCoeffs<T>::get(ellipse.getCore(), 0.0);
    cImage = afw::math::offsetImage(
        *cImage,
        ellipse.getCenter().getX(),
        ellipse.getCenter().getY(),
        ctrl.shiftKernel
    );
    if (!image.getBBox(afw::image::PARENT).contains(cImage->getBBox(afw::image::PARENT))) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % cImage->getBBox(afw::image::PARENT)).str()
        );
    }
    afw::image::Image<T> subImage(image, cImage->getBBox(afw::image::PARENT), afw::image::PARENT);
    return (subImage.getArray().template asEigen<Eigen::ArrayXpr>()
            * cImage->getArray().template asEigen<Eigen::ArrayXpr>()).sum();
}
template <typename T>
FluxComponent ApertureFluxAlgorithm::computeSincFlux(
    afw::image::MaskedImage<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    CONST_PTR(afw::image::Image<T>) cImage = SincCoeffs<T>::get(ellipse.getCore(), 0.0);
    cImage = afw::math::offsetImage(
        *cImage,
        ellipse.getCenter().getX(),
        ellipse.getCenter().getY(),
        ctrl.shiftKernel
    );
    if (!image.getBBox(afw::image::PARENT).contains(cImage->getBBox(afw::image::PARENT))) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % cImage->getBBox(afw::image::PARENT)).str()
        );
    }
    afw::image::MaskedImage<T> subImage(image, cImage->getBBox(afw::image::PARENT), afw::image::PARENT);
    FluxComponent result;
    result.flux = (subImage.getImage()->getArray().template asEigen<Eigen::ArrayXpr>()
                   * cImage->getArray().template asEigen<Eigen::ArrayXpr>()).sum();
    result.fluxSigma = std::sqrt(
        (subImage.getVariance()->getArray().template asEigen<Eigen::ArrayXpr>().template cast<T>()
         * cImage->getArray().template asEigen<Eigen::ArrayXpr>().square()).sum()
    );
    return result;
}

template <typename T>
Flux ApertureFluxAlgorithm::computeNaiveFlux(
    afw::image::Image<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    afw::geom::ellipses::PixelRegion region(ellipse); // behaves mostly like a Footprint
    if (!image.getBBox(afw::image::PARENT).contains(region.getBBox())) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % region.getBBox()).str()
        );
    }
    Flux flux = 0;
    for (
        afw::geom::ellipses::PixelRegion::Iterator spanIter = region.begin(), spanEnd = region.end();
        spanIter != spanEnd;
        ++spanIter
    ) {
        typename afw::image::Image<T>::x_iterator pixIter = image.x_at(
            spanIter->getBeginX() - image.getX0(),
            spanIter->getY() - image.getY0()
        );
        flux += std::accumulate(pixIter, pixIter + spanIter->getWidth(), 0.0);
    }
    return flux;
}

template <typename T>
FluxComponent ApertureFluxAlgorithm::computeNaiveFlux(
    afw::image::MaskedImage<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    FluxComponent result;
    afw::geom::ellipses::PixelRegion region(ellipse); // behaves mostly like a Footprint
    if (!image.getBBox(afw::image::PARENT).contains(region.getBBox())) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % region.getBBox()).str()
        );
    }
    result.flux = 0.0;
    result.fluxSigma = 0.0;
    for (
        afw::geom::ellipses::PixelRegion::Iterator spanIter = region.begin(), spanEnd = region.end();
        spanIter != spanEnd;
        ++spanIter
    ) {
        typename afw::image::MaskedImage<T>::Image::x_iterator pixIter = image.getImage()->x_at(
            spanIter->getBeginX() - image.getX0(),
            spanIter->getY() - image.getY0()
        );
        typename afw::image::MaskedImage<T>::Variance::x_iterator varIter = image.getVariance()->x_at(
            spanIter->getBeginX() - image.getX0(),
            spanIter->getY() - image.getY0()
        );
        result.flux += std::accumulate(pixIter, pixIter + spanIter->getWidth(), 0.0);
        // we use this to hold variance as we accumulate...
        result.fluxSigma += std::accumulate(varIter, varIter + spanIter->getWidth(), 0.0);
    }
    result.fluxSigma = std::sqrt(result.fluxSigma); // ...and switch back to sigma here.
    return result;
}


#define INSTANTIATE(T)                                          \
    template                                                    \
    Flux ApertureFluxAlgorithm::computeSincFlux(                \
        afw::image::Image<T> const &,                           \
        afw::geom::ellipses::Ellipse const &,                   \
        Control const &                                         \
    );                                                          \
    template                                                    \
    FluxComponent ApertureFluxAlgorithm::computeSincFlux(       \
        afw::image::MaskedImage<T> const &,                     \
        afw::geom::ellipses::Ellipse const &,                   \
        Control const &                                         \
    );                                                          \
    template                                                    \
    Flux ApertureFluxAlgorithm::computeNaiveFlux(               \
        afw::image::Image<T> const &,                           \
        afw::geom::ellipses::Ellipse const &,                   \
        Control const &                                         \
    );                                                          \
    template                                                    \
    FluxComponent ApertureFluxAlgorithm::computeNaiveFlux(      \
        afw::image::MaskedImage<T> const &,                     \
        afw::geom::ellipses::Ellipse const &,                   \
        Control const &                                         \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base
