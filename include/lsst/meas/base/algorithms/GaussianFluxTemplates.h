// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "boost/make_shared.hpp"
#include "boost/tuple/tuple.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/base/algorithms/SdssShapeImpl.h"
//#include "lsst/meas/base/algorithms/ScaledFlux.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwCoord = lsst::afw::coord;

namespace lsst {
namespace meas {
namespace base {
namespace algorithms {


/************************************************************************************************************/
    
template<typename ImageT>
std::pair<double, double>
getGaussianFlux(
    ImageT const& mimage,           // the data to process
    double background,               // background level
    double xcen, double ycen,         // centre of object
    double shiftmax,                  // max allowed centroid shift
    int maxIter=SDSS_SHAPE_MAX_ITER, ///< Maximum number of iterations
    float tol1=SDSS_SHAPE_TOL1, ///< Convergence tolerance for e1,e2
    float tol2=SDSS_SHAPE_TOL2, ///< Convergence tolerance for FWHM
    bool debug = false,
    PTR(SdssShapeImpl) shape=PTR(SdssShapeImpl)() // SDSS shape measurement
) {
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!shape) {
        shape = boost::make_shared<SdssShapeImpl>();
    }

    if (!getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, shape.get(),
                                    maxIter, tol1, tol2, debug)) {
        if (debug) {
            std::cout << "Not getAdaptiveMoments\n";
            for (int n = 0; n < SdssShapeImpl::N_FLAGS; ++n) {
                std::cout << shape->getFlag(SdssShapeImpl::Flag(n)) << " ";
            }
            std::cout << "\n";
        }
    } else {
        if (debug) {
            std::cout << "Not getAdaptiveMoments\n";
            for (int n = 0; n < SdssShapeImpl::N_FLAGS; ++n) {
                std::cout << shape->getFlag(SdssShapeImpl::Flag(n)) << " ";
            }
            std::cout << "\n";
        }
        double const scale = shape->getFluxScale();
        flux = scale*shape->getI0();
        fluxErr = scale*shape->getI0Err();
        if (debug) {
            std::cout << flux << " from getAdaptiveMoments\n";
        }
    }

    return std::make_pair(flux, fluxErr);
}




/*
 * Apply the algorithm to the PSF model
 */
double getPsfFactor(afwDet::Psf const & psf, afw::geom::Point2D const & center, double shiftmax,
                    int maxIter=SDSS_SHAPE_MAX_ITER, float tol1=SDSS_SHAPE_TOL1,
                    float tol2=SDSS_SHAPE_TOL2) {

    typedef afwDet::Psf::Image PsfImageT;
    PTR(PsfImageT) psfImage; // the image of the PSF
    PTR(PsfImageT) psfImageNoPad;   // Unpadded image of PSF
    
    int const pad = 5;
    try {
        psfImageNoPad = psf.computeImage(center);
        
        psfImage = PTR(PsfImageT)(
            new PsfImageT(psfImageNoPad->getDimensions() + afwGeom::Extent2I(2*pad))
            );
        afwGeom::BoxI middleBBox(afwGeom::Point2I(pad, pad), psfImageNoPad->getDimensions());
        
        PTR(PsfImageT) middle(new PsfImageT(*psfImage, middleBBox, afwImage::LOCAL));
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

}}}} // end lsst::meas::base::algorithms namespace
