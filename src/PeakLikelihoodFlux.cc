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
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/meas/base/PsfAttributes.h"

#include "lsst/meas/base/PeakLikelihoodFlux.h"

namespace lsst { namespace meas { namespace base {

PeakLikelihoodFluxAlgorithm::ResultMapper PeakLikelihoodFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

/**
Compute the value of one pixel of an image after a fractional pixel shift
Since we only want the value at one pixel, there is no need to shift the entire image;
instead we simply convolve at one point.

@throw pex::exceptions::RangeError if abs(fracShift) > 1 in either dimension
*/
template<typename T>
typename afw::image::MaskedImage<T>::SinglePixel computeShiftedValue(
    afw::image::MaskedImage<T> const &maskedImage, ///< masked image
    std::string const &warpingKernelName,   ///< warping kernel name
    afw::geom::Point2D const &fracShift,    ///< amount of sub-pixel shift (pixels)
    afw::geom::Point2I const &parentInd     ///< parent index at which to compute pixel
) {
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename afw::image::Image<double> KernelImageT;

    PTR(afw::math::SeparableKernel) warpingKernelPtr = afw::math::makeWarpingKernel(warpingKernelName);

    if ((std::abs(fracShift[0]) >= 1) || (std::abs(fracShift[1]) >= 1)) {
        std::ostringstream os;
        os << "fracShift = " << fracShift << " too large; abs value must be < 1 in both axes";
        throw LSST_EXCEPT(pex::exceptions::RangeError, os.str());
    }

    // warping kernels have even dimension and want the peak to the right of center
    if (fracShift[0] < 0) {
        warpingKernelPtr->setCtrX(warpingKernelPtr->getCtrX() + 1);
    }
    if (fracShift[1] < 0) {
        warpingKernelPtr->setCtrY(warpingKernelPtr->getCtrY() + 1);
    }
    afw::geom::Box2I warpingOverlapBBox(
        parentInd - afw::geom::Extent2I(warpingKernelPtr->getCtr()),
        warpingKernelPtr->getDimensions());
    if (!maskedImage.getBBox(afw::image::PARENT).contains(warpingOverlapBBox)) {
        std::ostringstream os;
        os << "Warping kernel extends off the edge"
            << "; kernel bbox = " << warpingOverlapBBox
            << "; exposure bbox = " << maskedImage.getBBox(afw::image::PARENT);
        throw LSST_EXCEPT(pex::exceptions::RangeError, os.str());
    }
    warpingKernelPtr->setKernelParameters(std::make_pair(fracShift[0], fracShift[1]));
    KernelImageT warpingKernelImage(warpingKernelPtr->getDimensions());
    warpingKernelPtr->computeImage(warpingKernelImage, true);
    typename KernelImageT::const_xy_locator const warpingKernelLoc = warpingKernelImage.xy_at(0,0);

    // Compute imLoc: an image locator that matches kernel locator (0,0) such that
    // image ctrPix overlaps center of warping kernel
    afw::geom::Point2I subimMin = warpingOverlapBBox.getMin();
    typename MaskedImageT::const_xy_locator const mimageLoc = maskedImage.xy_at(subimMin.getX(), subimMin.getY());
    return afw::math::convolveAtAPoint<MaskedImageT, MaskedImageT>(
        mimageLoc, warpingKernelLoc, warpingKernelPtr->getWidth(), warpingKernelPtr->getHeight());
}

template <typename T>
void PeakLikelihoodFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Result & result,
    Control const & ctrl
) {
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    MaskedImageT const& mimage = exposure.getMaskedImage();

    // BBox for data image
    afw::geom::BoxI imageBBox(mimage.getBBox(afw::image::PARENT));
/**
 * Given an image and a pixel position, return a Flux
 *
 * @throw lsst::pex::exceptions::InvalidParameterError if the exposure has no PSF.
 * @throw lsst::pex::exceptions::RangeError if the warping (centering) kernel
 *      is not fully contained within the exposure.
 * @throw lsst::pex::exceptions::RangeError if the center not within exposure.
 *      (This avoids insane center values from confusing the test for warping kernel within exposure).
 */

    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "exposure has no PSF");
    }
    PTR(afw::detection::Psf const) psfPtr = exposure.getPsf();
    if (!afw::geom::Box2D(mimage.getBBox(afw::image::PARENT)).contains(center)) {
        std::ostringstream os;
        os << "Center = " << center << " not in exposure bbox" << mimage.getBBox(afw::image::PARENT);
        throw LSST_EXCEPT(pex::exceptions::RangeError, os.str());
    }

    // compute parent index and fractional offset of ctrPix: the pixel closest to "center",
    // the centroid of the source
    std::pair<int, double> const xCtrPixParentIndFrac = afw::image::positionToIndex(center.getX(), true);
    std::pair<int, double> const yCtrPixParentIndFrac = afw::image::positionToIndex(center.getY(), true);

    afw::geom::Point2I ctrPixParentInd(xCtrPixParentIndFrac.first, yCtrPixParentIndFrac.first);
    afw::geom::Point2D ctrPixPos(
        afw::image::indexToPosition(ctrPixParentInd[0]),
        afw::image::indexToPosition(ctrPixParentInd[1])
    );

    // compute weight = 1/sum(PSF^2) for PSF at ctrPix, where PSF is normalized to a sum of 1
    PsfAttributes psfAttr(psfPtr, ctrPixParentInd);
    double weight = psfAttr.computeEffectiveArea();

    /*
     * Compute value of image at center of source, as shifted by a fractional pixel to center the source
     * on ctrPix.
     */
    typename MaskedImageT::SinglePixel mimageCtrPix = computeShiftedValue(
        mimage,
        ctrl.warpingKernelName,
        afw::geom::Point2D(xCtrPixParentIndFrac.second, yCtrPixParentIndFrac.second),
        ctrPixParentInd
    );
    double flux = mimageCtrPix.image()*weight;
    double var = mimageCtrPix.variance()*weight*weight;
    result.flux = flux;
    result.fluxSigma = std::sqrt(var);

    /* ******************************************************* */
}

template <typename T>
void PeakLikelihoodFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure, inputs.position, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void PeakLikelihoodFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void PeakLikelihoodFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base
