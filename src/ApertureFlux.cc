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
#include "lsst/afw/table/Source.h"
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

ApertureFluxAlgorithm::ApertureFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxKey(
        afw::table::ArrayKey<Flux>::addFields(
            schema,
            name + "_flux",
            "flux within %f pixel aperture",
            "dn",
            ctrl.radii
        )
    ),
    _fluxSigmaKey(
        afw::table::ArrayKey<Flux>::addFields(
            schema,
            name + "_fluxSigma",
            "1-sigma uncertainty on flux within %f pixel aperture",
            "dn",
            ctrl.radii
        )
    )
{
}

ApertureFluxAlgorithm::FlagKeys::FlagKeys(
    std::string const & name, afw::table::Schema & schema, int index
) :
    failed(
        schema.addField<afw::table::Flag>(
            (boost::format("%s_flag_%d") % name % index).str(),
            (boost::format("flag set if aperture %d failed for any reason") % index).str()
        )
    ),
    apertureTruncated(
        schema.addField<afw::table::Flag>(
            (boost::format("%s_flag_apertureTruncated_%d") % name % index).str(),
            (boost::format("flag set if aperture %d did not fit within the measurement image") % index).str()
        )
    ),
    sincCoeffsTruncated(
        schema.addField<afw::table::Flag>(
            (boost::format("%s_flag_sincCoeffsTruncated_%d") % name % index).str(),
            (boost::format("flag set if the full sinc coefficient image for aperture %d did not "
                           "fit within the measurement image") % index).str()
        )
    )
{}

void ApertureFluxAlgorithm::copyResultToRecord(
    Result const & result,
    afw::table::SourceRecord & record,
    int index
) const {
    record.set(_fluxKey[index], result.flux);
    record.set(_fluxSigmaKey[index], result.fluxSigma);
    if (result.getFlag(APERTURE_TRUNCATED)) {
        record.set(_flagKeys[index].apertureTruncated, true);
        record.set(_flagKeys[index].failed, true);
    }
    if (result.getFlag(SINC_COEFFS_TRUNCATED)) {
        record.set(_flagKeys[index].sincCoeffsTruncated, true);
        // TODO DM-464: set suspect flag
    }
}

namespace {

// Helper function for computeSincFlux get Sinc flux coefficients, and handle cases where the coeff
// image needs to be clipped to fit in the measurement image
template <typename T>
CONST_PTR(afw::image::Image<T>) getSincCoeffs(
    afw::geom::Box2I const & bbox,                // measurement image bbox we need to fit inside
    afw::geom::ellipses::Ellipse const & ellipse, // ellipse that defines the aperture
    ApertureFluxAlgorithm::Result & result,       // result object where we set flags if we do clip
    ApertureFluxAlgorithm::Control const & ctrl   // configuration
) {
    CONST_PTR(afw::image::Image<T>) cImage = SincCoeffs<T>::get(ellipse.getCore(), 0.0);
    cImage = afw::math::offsetImage(
        *cImage,
        ellipse.getCenter().getX(),
        ellipse.getCenter().getY(),
        ctrl.shiftKernel
    );
    if (!bbox.contains(cImage->getBBox())) {
        // We had to clip out at least part part of the coeff image,
        // but since that's much larger than the aperture (and close
        // to zero outside the aperture), it may not be a serious
        // problem.
        result.setFlag(ApertureFluxAlgorithm::SINC_COEFFS_TRUNCATED);
        afw::geom::Box2I overlap = cImage->getBBox();
        overlap.clip(bbox);
        if (!overlap.contains(afw::geom::Box2I(ellipse.computeBBox()))) {
            // The clipping was indeed serious, as we we did have to clip within
            // the aperture; can't expect any decent answer at this point.
            result.setFlag(ApertureFluxAlgorithm::APERTURE_TRUNCATED);
        }
        cImage = boost::make_shared< afw::image::Image<T> >(*cImage, overlap);
    }
    return cImage;
}

} // anonymous

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeSincFlux(
    afw::image::Image<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    Result result;
    CONST_PTR(afw::image::Image<T>) cImage = getSincCoeffs<T>(image.getBBox(), ellipse, result, ctrl);
    if (result.getFlag(APERTURE_TRUNCATED)) return result;
    afw::image::Image<T> subImage(image, cImage->getBBox());
    result.flux = (subImage.getArray().template asEigen<Eigen::ArrayXpr>()
                   * cImage->getArray().template asEigen<Eigen::ArrayXpr>()).sum();
    return result;
}

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeSincFlux(
    afw::image::MaskedImage<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    Result result;
    CONST_PTR(afw::image::Image<T>) cImage = getSincCoeffs<T>(image.getBBox(), ellipse, result, ctrl);
    if (result.getFlag(APERTURE_TRUNCATED)) return result;
    afw::image::MaskedImage<T> subImage(image, cImage->getBBox(afw::image::PARENT), afw::image::PARENT);
    result.flux = (subImage.getImage()->getArray().template asEigen<Eigen::ArrayXpr>()
                   * cImage->getArray().template asEigen<Eigen::ArrayXpr>()).sum();
    result.fluxSigma = std::sqrt(
        (subImage.getVariance()->getArray().template asEigen<Eigen::ArrayXpr>().template cast<T>()
         * cImage->getArray().template asEigen<Eigen::ArrayXpr>().square()).sum()
    );
    return result;
}

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeNaiveFlux(
    afw::image::Image<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    Result result;
    afw::geom::ellipses::PixelRegion region(ellipse); // behaves mostly like a Footprint
    if (!image.getBBox().contains(region.getBBox())) {
        result.setFlag(APERTURE_TRUNCATED);
        return result;
    }
    result.flux = 0;
    for (
        afw::geom::ellipses::PixelRegion::Iterator spanIter = region.begin(), spanEnd = region.end();
        spanIter != spanEnd;
        ++spanIter
    ) {
        typename afw::image::Image<T>::x_iterator pixIter = image.x_at(
            spanIter->getBeginX() - image.getX0(),
            spanIter->getY() - image.getY0()
        );
        result.flux += std::accumulate(pixIter, pixIter + spanIter->getWidth(), 0.0);
    }
    return result;
}

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeNaiveFlux(
    afw::image::MaskedImage<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    Result result;
    afw::geom::ellipses::PixelRegion region(ellipse); // behaves mostly like a Footprint
    if (!image.getBBox().contains(region.getBBox())) {
        result.setFlag(APERTURE_TRUNCATED);
        return result;
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

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeFlux(
    afw::image::Image<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    return (afw::geom::ellipses::Axes(ellipse.getCore()).getB() <= ctrl.maxSincRadius)
        ? computeSincFlux(image, ellipse, ctrl)
        : computeNaiveFlux(image, ellipse, ctrl);
}

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeFlux(
    afw::image::MaskedImage<T> const & image,
    afw::geom::ellipses::Ellipse const & ellipse,
    Control const & ctrl
) {
    return (afw::geom::ellipses::Axes(ellipse.getCore()).getB() <= ctrl.maxSincRadius)
        ? computeSincFlux(image, ellipse, ctrl)
        : computeNaiveFlux(image, ellipse, ctrl);
}
#define INSTANTIATE(T)                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeFlux( \
        afw::image::Image<T> const &,                                   \
        afw::geom::ellipses::Ellipse const &,                           \
        Control const &                                                 \
    );                                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeFlux( \
        afw::image::MaskedImage<T> const &,                             \
        afw::geom::ellipses::Ellipse const &,                           \
        Control const &                                                 \
    );                                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeSincFlux( \
        afw::image::Image<T> const &,                                   \
        afw::geom::ellipses::Ellipse const &,                           \
        Control const &                                                 \
    );                                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeSincFlux( \
        afw::image::MaskedImage<T> const &,                             \
        afw::geom::ellipses::Ellipse const &,                           \
        Control const &                                                 \
    );                                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeNaiveFlux( \
        afw::image::Image<T> const &,                                   \
        afw::geom::ellipses::Ellipse const &,                           \
        Control const &                                                 \
    );                                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::computeNaiveFlux( \
        afw::image::MaskedImage<T> const &,                             \
        afw::geom::ellipses::Ellipse const &,                           \
        Control const &                                                 \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base
src/CircularApertureFlux.cc                                                                         0000644 0001756 0000144 00000005231 12434274776 015155  0                                                                                                    ustar   pgee                            users                                                                                                                                                                                                                  // -*- lsst-c++ -*-
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

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/CircularApertureFlux.h"
#include "lsst/meas/base/SincCoeffs.h"

namespace lsst { namespace meas { namespace base {

CircularApertureFluxAlgorithm::CircularApertureFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema,
    daf::base::PropertySet & metadata
) : ApertureFluxAlgorithm(ctrl, name, schema),
    _centroidExtractor(schema, name)
{
    for (std::size_t i = 0; i < ctrl.radii.size(); ++i) {
        if (ctrl.radii[i] > ctrl.maxSincRadius) break;
        SincCoeffs<float>::cache(0.0, ctrl.radii[i]);
        metadata.add(name + "_radii", ctrl.radii[i]);
    }
    _flagKeys.reserve(ctrl.radii.size());
    for (std::size_t i = 0; i < ctrl.radii.size(); ++i) {
        _flagKeys.push_back(FlagKeys(name, schema, i));
    }
}

void CircularApertureFluxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    afw::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);   
    afw::geom::ellipses::Ellipse ellipse(afw::geom::ellipses::Axes(1.0, 1.0, 0.0), center);
    PTR(afw::geom::ellipses::Axes) axes
        = boost::static_pointer_cast<afw::geom::ellipses::Axes>(ellipse.getCorePtr());
    for (std::size_t i = 0; i < _ctrl.radii.size(); ++i) {
        axes->setA(_ctrl.radii[i]);
        axes->setB(_ctrl.radii[i]);
        ApertureFluxAlgorithm::Result result = computeFlux(exposure.getMaskedImage(), ellipse, _ctrl);
        copyResultToRecord(result, measRecord, i);
    }
}

void CircularApertureFluxAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}} // namespace lsst::meas::base
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       