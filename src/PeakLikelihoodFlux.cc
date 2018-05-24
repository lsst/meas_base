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
#include "lsst/geom/Box.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

#include "lsst/meas/base/PeakLikelihoodFlux.h"
#include "lsst/afw/table/Source.h"

namespace lsst {
namespace meas {
namespace base {
namespace {
FlagDefinitionList flagDefinitions;
}  // namespace

FlagDefinition const PeakLikelihoodFluxAlgorithm::FAILURE = flagDefinitions.addFailureFlag();

FlagDefinitionList const &PeakLikelihoodFluxAlgorithm::getFlagDefinitions() { return flagDefinitions; }

namespace {

/************************************************************************************************************/
/**
 * @class PsfAttributes
 *
 * A class to contain various attributes of the Psf
 * - most notably, a width (1-D RMS size) to be used to
 *   make a single gaussian psf for fast convolution.
 *
 * \deprecated
 * This class is deprecated in favour of virtual methods on Psf
 *
 * An example of the new API is:
 * \code
 * afw::geom::ellipses::Quadrupole shape = psf->computeShape();
 * double const smoothingSigma = shape.getDeterminantRadius();
 * \endcode
 */
class PsfAttributes {
public:
    enum Method {
        ADAPTIVE_MOMENT,   ///< Calculate width using adaptive Gaussian weights
        FIRST_MOMENT,      ///< Calculate width using \<r>
        SECOND_MOMENT,     ///< Calculate width using \<r^2>
        NOISE_EQUIVALENT,  ///< Calculate width as sqrt(n_eff/(4 pi))
        BICKERTON          ///< Weight \<r^2> by I^2 to avoid negative fluxes
    };

    PsfAttributes(CONST_PTR(afw::detection::Psf) psf, int const iX, int const iY);
    PsfAttributes(CONST_PTR(afw::detection::Psf) psf, geom::Point2I const &cen);

    double computeGaussianWidth(Method how = ADAPTIVE_MOMENT) const;
    double computeEffectiveArea() const;

private:
    PTR(afw::image::Image<double>) _psfImage;
};

/**
 * @brief Constructor for PsfAttributes
 */
PsfAttributes::PsfAttributes(CONST_PTR(afw::detection::Psf) psf,  ///< The psf whose attributes we want
                             int const iX,  ///< the x position in the frame we want the attributes at
                             int const iY   ///< the y position in the frame we want the attributes at
) {
    // N.b. (iX, iY) are ints so that we know this image is centered in the central pixel of _psfImage
    _psfImage = psf->computeImage(geom::PointD(iX, iY));
}

/**
 * @brief Constructor for PsfAttributes
 */
PsfAttributes::PsfAttributes(
        CONST_PTR(afw::detection::Psf) psf,  ///< The psf whose attributes we want
        geom::Point2I const &cen             ///< the position in the frame we want the attributes at
        )
        :  // N.b. cen is a PointI so that we know this image is centered in the central pixel of _psfImage
          _psfImage(psf->computeImage(geom::PointD(cen))) {}

/**
 * @brief Compute the effective area of the psf ( sum(I)^2/sum(I^2) )
 *
 */
double PsfAttributes::computeEffectiveArea() const {
    double sum = 0.0;
    double sumsqr = 0.0;
    for (int iY = 0; iY != _psfImage->getHeight(); ++iY) {
        afw::image::Image<double>::x_iterator end = _psfImage->row_end(iY);
        for (afw::image::Image<double>::x_iterator ptr = _psfImage->row_begin(iY); ptr != end; ++ptr) {
            sum += *ptr;
            sumsqr += (*ptr) * (*ptr);
        }
    }
    return sum * sum / sumsqr;
}

}  // end anonymous namespace

/**
Compute the value of one pixel of an image after a fractional pixel shift
Since we only want the value at one pixel, there is no need to shift the entire image;
instead we simply convolve at one point.

@throw pex::exceptions::RangeError if abs(fracShift) > 1 in either dimension
*/
template <typename T>
typename afw::image::MaskedImage<T>::SinglePixel computeShiftedValue(
        afw::image::MaskedImage<T> const &maskedImage,  ///< masked image
        std::string const &warpingKernelName,           ///< warping kernel name
        geom::Point2D const &fracShift,                 ///< amount of sub-pixel shift (pixels)
        geom::Point2I const &parentInd                  ///< parent index at which to compute pixel
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
    geom::Box2I warpingOverlapBBox(parentInd - geom::Extent2I(warpingKernelPtr->getCtr()),
                                   warpingKernelPtr->getDimensions());
    if (!maskedImage.getBBox().contains(warpingOverlapBBox)) {
        std::ostringstream os;
        os << "Warping kernel extends off the edge"
           << "; kernel bbox = " << warpingOverlapBBox << "; exposure bbox = " << maskedImage.getBBox();
        throw LSST_EXCEPT(pex::exceptions::RangeError, os.str());
    }
    warpingKernelPtr->setKernelParameters(std::make_pair(fracShift[0], fracShift[1]));
    KernelImageT warpingKernelImage(warpingKernelPtr->getDimensions());
    warpingKernelPtr->computeImage(warpingKernelImage, true);
    typename KernelImageT::const_xy_locator const warpingKernelLoc = warpingKernelImage.xy_at(0, 0);

    // Compute imLoc: an image locator that matches kernel locator (0,0) such that
    // image ctrPix overlaps center of warping kernel
    geom::Point2I subimMin = warpingOverlapBBox.getMin();
    typename MaskedImageT::const_xy_locator const mimageLoc =
            maskedImage.xy_at(subimMin.getX(), subimMin.getY());
    return afw::math::convolveAtAPoint<MaskedImageT, MaskedImageT>(
            mimageLoc, warpingKernelLoc, warpingKernelPtr->getWidth(), warpingKernelPtr->getHeight());
}
PeakLikelihoodFluxAlgorithm::PeakLikelihoodFluxAlgorithm(Control const &ctrl, std::string const &name,
                                                         afw::table::Schema &schema)
        : _ctrl(ctrl),
          _fluxResultKey(FluxResultKey::addFields(schema, name, "flux from PeakLikelihood Flux algorithm")),
          _centroidExtractor(schema, name) {
    _flagHandler = FlagHandler::addFields(schema, name, getFlagDefinitions());
}

void PeakLikelihoodFluxAlgorithm::measure(afw::table::SourceRecord &measRecord,
                                          afw::image::Exposure<float> const &exposure) const {
    // get the value from the centroid slot only
    geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    FluxResult result;
    typedef afw::image::Exposure<float>::MaskedImageT MaskedImageT;
    MaskedImageT const &mimage = exposure.getMaskedImage();

    /**
     * Given an image and a pixel position, return a Flux
     *
     * @throw pex::exceptions::InvalidParameterError if the exposure has no PSF.
     * @throw pex::exceptions::RangeError if the warping (centering) kernel
     *      is not fully contained within the exposure.
     * @throw pex::exceptions::RangeError if the center not within exposure.
     *      (This avoids insane center values from confusing the test for warping kernel within exposure).
     */

    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "exposure has no PSF");
    }
    PTR(afw::detection::Psf const) psfPtr = exposure.getPsf();
    if (!geom::Box2D(mimage.getBBox()).contains(center)) {
        std::ostringstream os;
        os << "Center = " << center << " not in exposure bbox" << mimage.getBBox();
        throw LSST_EXCEPT(pex::exceptions::RangeError, os.str());
    }

    // compute parent index and fractional offset of ctrPix: the pixel closest to "center",
    // the centroid of the source
    std::pair<int, double> const xCtrPixParentIndFrac = afw::image::positionToIndex(center.getX(), true);
    std::pair<int, double> const yCtrPixParentIndFrac = afw::image::positionToIndex(center.getY(), true);

    geom::Point2I ctrPixParentInd(xCtrPixParentIndFrac.first, yCtrPixParentIndFrac.first);
    geom::Point2D ctrPixPos(afw::image::indexToPosition(ctrPixParentInd[0]),
                            afw::image::indexToPosition(ctrPixParentInd[1]));

    // compute weight = 1/sum(PSF^2) for PSF at ctrPix, where PSF is normalized to a sum of 1
    PsfAttributes psfAttr(psfPtr, ctrPixParentInd);
    double weight = psfAttr.computeEffectiveArea();

    /*
     * Compute value of image at center of source, as shifted by a fractional pixel to center the source
     * on ctrPix.
     */
    MaskedImageT::SinglePixel mimageCtrPix = computeShiftedValue(
            mimage, _ctrl.warpingKernelName,
            geom::Point2D(xCtrPixParentIndFrac.second, yCtrPixParentIndFrac.second), ctrPixParentInd);
    double flux = mimageCtrPix.image() * weight;
    double var = mimageCtrPix.variance() * weight * weight;
    result.flux = flux;
    result.fluxSigma = std::sqrt(var);
    measRecord.set(_fluxResultKey, result);
    _flagHandler.setValue(measRecord, FAILURE.number, false);
}

void PeakLikelihoodFluxAlgorithm::fail(afw::table::SourceRecord &measRecord, MeasurementError *error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
