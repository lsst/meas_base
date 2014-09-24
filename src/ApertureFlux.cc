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

ApertureFluxControl::ApertureFluxControl() : radii(10), maxSincRadius(10.0), shiftKernel("lanczos5") {
    // defaults here stolen from HSC pipeline defaults
    static boost::array<double,10> defaultRadii = {{
        3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0
    }};
    std::copy(defaultRadii.begin(), defaultRadii.end(), radii.begin());
}

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

ApertureFluxAlgorithm::ResultMapper ApertureFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema,
    std::string const & name,
    Control const & ctrl
) {
    return ResultMapper(schema, name, ctrl.radii);
}

template <typename T>
void ApertureFluxAlgorithm::apply(
    afw::image::MaskedImage<T> const & image,
    afw::geom::Point2D const & position,
    Result & result,
    Control const & ctrl
) {
    ApertureFluxComponent arrays(ctrl.radii.size()); // invoke this constructor to initialize arrays
    result.flux = arrays.flux;
    result.fluxSigma = arrays.fluxSigma;
    for (std::size_t n = 0; n < ctrl.radii.size(); ++n) {
        afw::geom::ellipses::Ellipse ellipse(
            afw::geom::ellipses::Axes(ctrl.radii[n], ctrl.radii[n]),
            position
        );
        result.set(n, computeFlux(image, ellipse, ctrl));
    }
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
            MeasurementError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % cImage->getBBox(afw::image::PARENT)).str(),
            EDGE
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
            MeasurementError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % cImage->getBBox(afw::image::PARENT)).str(),
            EDGE
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
            MeasurementError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % region.getBBox()).str(),
            EDGE
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
            MeasurementError,
            (boost::format("Measurement image with (bbox=%s) is not large enough for aperture (bbox=%s)")
             % image.getBBox(afw::image::PARENT) % region.getBBox()).str(),
            EDGE
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
    void ApertureFluxAlgorithm::apply(                          \
        afw::image::MaskedImage<T> const &,                     \
        afw::geom::Point2D const &,                             \
        Result &,                                               \
        Control const &                                         \
    );                                                          \
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
