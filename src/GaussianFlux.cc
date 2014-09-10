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

#include "Eigen/Core"
#include "Eigen/LU"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/meas/base/GaussianFlux.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {


/************************************************************************************************************/
    
template<typename ImageT>
std::pair<double, double>
getGaussianFlux(
    ImageT const& mimage,           // the data to process
    double background,               // background level
    double xcen, double ycen,         // centre of object
    double shiftmax,                  // max allowed centroid shift
    int maxIter=detail::SDSS_SHAPE_MAX_ITER, ///< Maximum number of iterations
    float tol1=detail::SDSS_SHAPE_TOL1, ///< Convergence tolerance for e1,e2
    float tol2=detail::SDSS_SHAPE_TOL2, ///< Convergence tolerance for FWHM
    PTR(detail::SdssShapeImpl) shape=PTR(detail::SdssShapeImpl)() // detail::SDSS shape measurement
) {
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!shape) {
        shape = boost::make_shared<detail::SdssShapeImpl>();
    }

    if (!getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, shape.get(),
                                    maxIter, tol1, tol2)) {
    } else {
        double const scale = shape->getFluxScale();
        flux = scale*shape->getI0();
        fluxErr = scale*shape->getI0Err();
    }

    return std::make_pair(flux, fluxErr);
}




/*
 * Apply the algorithm to the PSF model
 */
double getPsfFactor(lsst::afw::detection::Psf const & psf, afw::geom::Point2D const & center, double shiftmax,
                    int maxIter=detail::SDSS_SHAPE_MAX_ITER, float tol1=detail::SDSS_SHAPE_TOL1,
                    float tol2=detail::SDSS_SHAPE_TOL2) {

    typedef lsst::afw::detection::Psf::Image PsfImageT;
    PTR(PsfImageT) psfImage; // the image of the PSF
    PTR(PsfImageT) psfImageNoPad;   // Unpadded image of PSF
    
    int const pad = 5;
    try {
        psfImageNoPad = psf.computeImage(center);
        
        psfImage = PTR(PsfImageT)(
            new PsfImageT(psfImageNoPad->getDimensions() + lsst::afw::geom::Extent2I(2*pad))
            );
        lsst::afw::geom::BoxI middleBBox(lsst::afw::geom::Point2I(pad, pad), psfImageNoPad->getDimensions());
        
        PTR(PsfImageT) middle(new PsfImageT(*psfImage, middleBBox, lsst::afw::image::LOCAL));
        *middle <<= *psfImageNoPad;
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)")
                            % center.getX() % center.getY()).str());
        throw e;
    }
    // Estimate the GaussianFlux for the Psf
    double const psfXCen = 0.5*(psfImage->getWidth() - 1); // Center of (21x21) image is (10.0, 10.0)
    double const psfYCen = 0.5*(psfImage->getHeight() - 1);
    std::pair<double, double> const result = getGaussianFlux(*psfImage, 0.0, psfXCen, psfYCen, shiftmax,
                                                             maxIter, tol1, tol2);
    return result.first;
}


GaussianFluxAlgorithm::ResultMapper GaussianFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
void GaussianFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Result & result,
    Control const & ctrl
) {
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            getFlagDefinitions()[NO_PSF].doc,
            NO_PSF
        );
    }
    PTR(afw::detection::Psf::Image) psfImage = psf->computeImage(center);
    afw::geom::Box2I fitBBox = psfImage->getBBox(afw::image::PARENT);
    fitBBox.clip(exposure.getBBox(afw::image::PARENT));
    if (fitBBox != psfImage->getBBox(afw::image::PARENT)) {
        result.setFlag(EDGE);
    }
    afw::detection::Footprint fitRegion(fitBBox);
    if (!ctrl.badMaskPlanes.empty()) {
        afw::image::MaskPixel badBits = 0x0;
        for (
            std::vector<std::string>::const_iterator i = ctrl.badMaskPlanes.begin();
            i != ctrl.badMaskPlanes.end();
            ++i
        ) {
            badBits |= exposure.getMaskedImage().getMask()->getPlaneBitMask(*i);
        }
        fitRegion.intersectMask(*exposure.getMaskedImage().getMask(), badBits);
    }
    if (fitRegion.getArea() == 0) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            getFlagDefinitions()[NO_GOOD_PIXELS].doc,
            NO_GOOD_PIXELS
        );
    }
    //  This code came straight out of the GaussianFlux.apply() in meas_algorithms with few changes
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    typename afw::image::Exposure<T>::MaskedImageT const& mimage = exposure.getMaskedImage();

    double const xcen = center.getX() - mimage.getX0(); ///< column position in image pixel coords
    double const ycen = center.getY() - mimage.getY0(); ///< row position

    std::pair<double, double> oldresult;
    if (ctrl.fixed) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            getFlagDefinitions()[NO_FIXED].doc,
            NO_FIXED
        );
    /*
        // Fixed aperture, defined by detail::SDSS shape measurement made elsewhere
        if (source.get(_shapeFlagKey)) {
            throw LSST_EXCEPT(pexExceptions::RuntimeError, "Shape measurement failed");
        }
        detail::SdssShapeImpl sdss(source.get(_centroidKey), source.get(_shapeKey));
        oldresult = detail::getFixedMomentsFlux(mimage, ctrl.background, xcen, ycen, sdss);
    */
    } else {
        // FIXME: propagate detail::SDSS shape measurement flags.
        /*
         * Find the object's adaptive-moments.  N.b. it would be better to use the SdssShape measurement
         * as this code repeats the work of that measurement
         */
        oldresult = getGaussianFlux<MaskedImageT>(mimage, ctrl.background, xcen, ycen, ctrl.shiftmax, ctrl.maxIter,
                                 ctrl.tol1, ctrl.tol2);
    }

    result.flux =  oldresult.first;
    result.fluxSigma = oldresult.second;

/*  Remove the psf scaling for first port -- pgee
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            lsst::meas::base::MeasurementError,
            getFlagDefinitions()[NO_PSF].doc,
            NO_PSF
        );
    }
    double psfFactor = getPsfFactor(*exposure.getPsf(), center, ctrl.shiftmax,
                                    ctrl.maxIter, ctrl.tol1, ctrl.tol2);
    result.fluxCorrectionKeys.psfFactor =  psfFactor;
*/
    //  End of meas_algorithms code
}

template <typename T>
void GaussianFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure, inputs.position, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void GaussianFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void GaussianFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

