// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include <cmath>
#include <tuple>

#include "boost/tuple/tuple.hpp"
#include "Eigen/LU"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/table/Source.h"

#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/base/SdssShape.h"

namespace pexPolicy = lsst::pex::policy;
namespace pexExceptions = lsst::pex::exceptions;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwTable = lsst::afw::table;

namespace lsst { namespace meas { namespace base {
namespace {
FlagDefinitionList flagDefinitions;
} // end anonymous

FlagDefinition const SdssShapeAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
FlagDefinition const SdssShapeAlgorithm::UNWEIGHTED_BAD = flagDefinitions.add("flag_unweightedBad", "Both weighted and unweighted moments were invalid");
FlagDefinition const SdssShapeAlgorithm::UNWEIGHTED = flagDefinitions.add("flag_unweighted", "Weighted moments converged to an invalid value; using unweighted moments");
FlagDefinition const SdssShapeAlgorithm::SHIFT = flagDefinitions.add("flag_shift", "centroid shifted by more than the maximum allowed amount");
FlagDefinition const SdssShapeAlgorithm::MAXITER = flagDefinitions.add("flag_maxIter", "Too many iterations in adaptive moments");
FlagDefinition const SdssShapeAlgorithm::PSF_SHAPE_BAD = flagDefinitions.add("flag_psf", "Failure in measuring PSF model shape");

FlagDefinitionList const & SdssShapeAlgorithm::getFlagDefinitions() {
    return flagDefinitions;
}


namespace {  // anonymous

typedef Eigen::Matrix<double,4,4,Eigen::DontAlign> Matrix4d;

// Return multiplier that transforms I0 to a total flux
double computeFluxScale(SdssShapeResult const & result) {
    /*
     * The shape is an ellipse that's axis-aligned in (u, v) [<uv> = 0] after rotation by theta:
     * <x^2> + <y^2> = <u^2> + <v^2>
     * <x^2> - <y^2> = cos(2 theta)*(<u^2> - <v^2>)
     * 2*<xy>        = sin(2 theta)*(<u^2> - <v^2>)
     */
    double const Mxx = result.xx; // <x^2>
    double const Mxy = result.xy; // <xy>
    double const Myy = result.yy; // <y^2>

    double const Muu_p_Mvv = Mxx + Myy;                             // <u^2> + <v^2>
    double const Muu_m_Mvv = ::sqrt(::pow(Mxx - Myy, 2) + 4*::pow(Mxy, 2)); // <u^2> - <v^2>
    double const Muu = 0.5*(Muu_p_Mvv + Muu_m_Mvv);
    double const Mvv = 0.5*(Muu_p_Mvv - Muu_m_Mvv);

    return lsst::afw::geom::TWOPI * ::sqrt(Muu*Mvv);
}

/*****************************************************************************/
/*
 * Error analysis, courtesy of David Johnston, University of Chicago
 */
/*
 * This function takes the 4 Gaussian parameters A, sigmaXXW and the
 * sky variance and fills in the Fisher matrix from the least squares fit.
 *
 * Following "Numerical Recipes in C" section 15.5, it ignores the 2nd
 * derivative parts and so the fisher matrix is just a function of these
 * best fit model parameters. The components are calculated analytically.
 */
Matrix4d
calc_fisher(SdssShapeResult const& shape, // the Shape that we want the the Fisher matrix for
            float bkgd_var              // background variance level for object
) {
    float const A = shape.flux;     // amplitude; will be converted to flux later
    float const sigma11W = shape.xx;
    float const sigma12W = shape.xy;
    float const sigma22W = shape.yy;

    double const D = sigma11W*sigma22W - sigma12W*sigma12W;

    if (D <= std::numeric_limits<double>::epsilon()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainError,
                          "Determinant is too small calculating Fisher matrix");
    }
/*
 * a normalization factor
 */
    if (bkgd_var <= 0.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainError,
                          (boost::format("Background variance must be positive (saw %g)") % bkgd_var).str());
    }
    double const F = afwGeom::PI*sqrt(D)/bkgd_var;
/*
 * Calculate the 10 independent elements of the 4x4 Fisher matrix
 */
    Matrix4d fisher;

    double fac = F*A/(4.0*D);
    fisher(0, 0) =  F;
    fisher(0, 1) =  fac*sigma22W;
    fisher(1, 0) =  fisher(0, 1);
    fisher(0, 2) =  fac*sigma11W;
    fisher(2, 0) =  fisher(0, 2);
    fisher(0, 3) = -fac*2*sigma12W;
    fisher(3, 0) =  fisher(0, 3);

    fac = 3.0*F*A*A/(16.0*D*D);
    fisher(1, 1) =  fac*sigma22W*sigma22W;
    fisher(2, 2) =  fac*sigma11W*sigma11W;
    fisher(3, 3) =  fac*4.0*(sigma12W*sigma12W + D/3.0);

    fisher(1, 2) =  fisher(3, 3)/4.0;
    fisher(2, 1) =  fisher(1, 2);
    fisher(1, 3) =  fac*(-2*sigma22W*sigma12W);
    fisher(3, 1) =  fisher(1, 3);
    fisher(2, 3) =  fac*(-2*sigma11W*sigma12W);
    fisher(3, 2) =  fisher(2, 3);

    return fisher;
}
//
// Here's a class to allow us to get the Image and variance from an Image or MaskedImage
//
template<typename ImageT>               // general case
struct ImageAdaptor {
    typedef ImageT Image;

    static bool const hasVariance = false;

    Image const& getImage(ImageT const& image) const {
        return image;
    }

    double getVariance(ImageT const&, int, int) {
        return std::numeric_limits<double>::quiet_NaN();
    }
};

template<typename T>                    // specialise to a MaskedImage
struct ImageAdaptor<afwImage::MaskedImage<T> > {
    typedef typename afwImage::MaskedImage<T>::Image Image;

    static bool const hasVariance = true;

    Image const& getImage(afwImage::MaskedImage<T> const& mimage) const {
        return *mimage.getImage();
    }

    double getVariance(afwImage::MaskedImage<T> const& mimage, int ix, int iy) {
        return mimage.at(ix, iy).variance();
    }
};

/// Calculate weights from moments
std::tuple<std::pair<bool, double>, double, double, double>
getWeights(double sigma11, double sigma12, double sigma22) {
    double const NaN = std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(sigma11) || std::isnan(sigma12) || std::isnan(sigma22)) {
        return std::make_tuple(std::make_pair(false, NaN), NaN, NaN, NaN);
    }
    double const det = sigma11*sigma22 - sigma12*sigma12; // determinant of sigmaXX matrix
    if (std::isnan(det) || det < std::numeric_limits<float>::epsilon()) { // a suitably small number
        return std::make_tuple(std::make_pair(false, det), NaN, NaN, NaN);
    }
    return std::make_tuple(std::make_pair(true, det), sigma22/det, -sigma12/det, sigma11/det);
}

/// Should we be interpolating?
bool shouldInterp(double sigma11, double sigma22, double det) {
    float const xinterp = 0.25; // I.e. 0.5*0.5
    return (sigma11 < xinterp || sigma22 < xinterp || det < xinterp*xinterp);
}

// Decide on the bounding box for the region to examine while calculating the adaptive moments
// This routine will work in either LOCAL or PARENT coordinates (but of course which you pass
// determine which you will get back).
afw::geom::Box2I computeAdaptiveMomentsBBox(
    afw::geom::Box2I const & bbox,  // full image bbox
    afw::geom::Point2D const & center, // centre of object
    double sigma11_w,              // quadratic moments of the
    double ,                       //         weighting function
    double sigma22_w,              //                    xx, xy, and yy
    double maxRadius = 1000        // Maximum radius of area to use
) {
    double radius = std::min(4*std::sqrt(std::max(sigma11_w, sigma22_w)), maxRadius);
    afw::geom::Extent2D offset(radius, radius);
    afw::geom::Box2I result(afw::geom::Box2D(center - offset, center + offset));
    result.clip(bbox);
    return result;
}

/*****************************************************************************/
/*
 * Calculate weighted moments of an object up to 2nd order
 */
template<bool fluxOnly, typename ImageT>
static int
calcmom(ImageT const& image,            // the image data
        float xcen, float ycen,         // centre of object
        lsst::afw::geom::BoxI bbox,    // bounding box to consider
        float bkgd,                     // data's background level
        bool interpflag,                // interpolate within pixels?
        double w11, double w12, double w22, // weights
        double *pI0,                        // amplitude of fit
        double *psum,                       // sum w*I (if !NULL)
        double *psumx, double *psumy,       // sum [xy]*w*I (if !fluxOnly)
        double *psumxx, double *psumxy, double *psumyy, // sum [xy]^2*w*I (if !fluxOnly)
        double *psums4,                                 // sum w*I*weight^2 (if !fluxOnly && !NULL)
        bool negative = false
       )
{

    float tmod, ymod;
    float X, Y;                          // sub-pixel interpolated [xy]
    float weight;
    float tmp;
    double sum, sumx, sumy, sumxx, sumyy, sumxy, sums4;
#define RECALC_W 0                      // estimate sigmaXX_w within BBox?
#if RECALC_W
    double wsum, wsumxx, wsumxy, wsumyy;

    wsum = wsumxx = wsumxy = wsumyy = 0;
#endif

    assert(w11 >= 0);                   // i.e. it was set
    if (fabs(w11) > 1e6 || fabs(w12) > 1e6 || fabs(w22) > 1e6) {
        return(-1);
    }

    sum = sumx = sumy = sumxx = sumxy = sumyy = sums4 = 0;

    int const ix0 = bbox.getMinX();       // corners of the box being analyzed
    int const ix1 = bbox.getMaxX();
    int const iy0 = bbox.getMinY();       // corners of the box being analyzed
    int const iy1 = bbox.getMaxY();

    if (ix0 < 0 || ix1 >= image.getWidth() || iy0 < 0 || iy1 >= image.getHeight()) {
        return -1;
    }

    for (int i = iy0; i <= iy1; ++i) {
        typename ImageT::x_iterator ptr = image.x_at(ix0, i);
        float const y = i - ycen;
        float const y2 = y*y;
        float const yl = y - 0.375;
        float const yh = y + 0.375;
        for (int j = ix0; j <= ix1; ++j, ++ptr) {
            float x = j - xcen;
            if (interpflag) {
                float const xl = x - 0.375;
                float const xh = x + 0.375;

                float expon = xl*xl*w11 + yl*yl*w22 + 2.0*xl*yl*w12;
                tmp = xh*xh*w11 + yh*yh*w22 + 2.0*xh*yh*w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xl*xl*w11 + yh*yh*w22 + 2.0*xl*yh*w12;
                expon = (expon > tmp) ? expon : tmp;
                tmp = xh*xh*w11 + yl*yl*w22 + 2.0*xh*yl*w12;
                expon = (expon > tmp) ? expon : tmp;

                if (expon <= 9.0) {
                    tmod = *ptr - bkgd;
                    for (Y = yl; Y <= yh; Y += 0.25) {
                        double const interpY2 = Y*Y;
                        for (X = xl; X <= xh; X += 0.25) {
                            double const interpX2 = X*X;
                            double const interpXy = X*Y;
                            expon = interpX2*w11 + 2*interpXy*w12 + interpY2*w22;
                            weight = std::exp(-0.5*expon);

                            ymod = tmod*weight;
                            sum += ymod;
                            if (!fluxOnly) {
                                sumx += ymod*(X + xcen);
                                sumy += ymod*(Y + ycen);
#if RECALC_W
                                wsum += weight;

                                tmp = interpX2*weight;
                                wsumxx += tmp;
                                sumxx += tmod*tmp;

                                tmp = interpXy*weight;
                                wsumxy += tmp;
                                sumxy += tmod*tmp;

                                tmp = interpY2*weight;
                                wsumyy += tmp;
                                sumyy += tmod*tmp;
#else
                                sumxx += interpX2*ymod;
                                sumxy += interpXy*ymod;
                                sumyy += interpY2*ymod;
#endif
                                sums4 += expon*expon*ymod;
                            }
                        }
                    }
                }
            } else {
                float x2 = x*x;
                float xy = x*y;
                float expon = x2*w11 + 2*xy*w12 + y2*w22;

                if (expon <= 14.0) {
                    weight = std::exp(-0.5*expon);
                    tmod = *ptr - bkgd;
                    ymod = tmod*weight;
                    sum += ymod;
                    if (!fluxOnly) {
                        sumx += ymod*j;
                        sumy += ymod*i;
#if RECALC_W
                        wsum += weight;

                        tmp = x2*weight;
                        wsumxx += tmp;
                        sumxx += tmod*tmp;

                        tmp = xy*weight;
                        wsumxy += tmp;
                        sumxy += tmod*tmp;

                        tmp = y2*weight;
                        wsumyy += tmp;
                        sumyy += tmod*tmp;
#else
                        sumxx += x2*ymod;
                        sumxy += xy*ymod;
                        sumyy += y2*ymod;
#endif
                        sums4 += expon*expon*ymod;
                    }
                }
            }
        }
    }


    std::tuple<std::pair<bool, double>, double, double, double> const weights = getWeights(w11, w12, w22);
    double const detW = std::get<1>(weights)*std::get<3>(weights) - std::pow(std::get<2>(weights), 2);
    *pI0 = sum/(afwGeom::PI*sqrt(detW));
    if (psum) {
        *psum = sum;
    }
    if (!fluxOnly) {
        *psumx = sumx;
        *psumy = sumy;
        *psumxx = sumxx;
        *psumxy = sumxy;
        *psumyy = sumyy;
        if (psums4 != NULL) {
            *psums4 = sums4;
        }
    }

#if RECALC_W
    if (wsum > 0 && !fluxOnly) {
        double det = w11*w22 - w12*w12;
        wsumxx /= wsum;
        wsumxy /= wsum;
        wsumyy /= wsum;
        printf("%g %g %g  %g %g %g\n", w22/det, -w12/det, w11/det, wsumxx, wsumxy, wsumyy);
    }
#endif

    if (negative) {
        return (fluxOnly || (sum < 0 && sumxx < 0 && sumyy < 0)) ? 0 : -1;
    } else {
        return (fluxOnly || (sum > 0 && sumxx > 0 && sumyy > 0)) ? 0 : -1;
    }
}

/*
 * Workhorse for adaptive moments
 *
 * All inputs are expected to be in LOCAL image coordinates
 */
template<typename ImageT>
bool getAdaptiveMoments(ImageT const& mimage, double bkgd, double xcen, double ycen, double shiftmax,
                        SdssShapeResult *shape, int maxIter, float tol1, float tol2, bool negative)
{
    double I0 = 0;                      // amplitude of best-fit Gaussian
    double sum;                         // sum of intensity*weight
    double sumx, sumy;                  // sum ((int)[xy])*intensity*weight
    double sumxx, sumxy, sumyy;         // sum {x^2,xy,y^2}*intensity*weight
    double sums4;                       // sum intensity*weight*exponent^2
    float const xcen0 = xcen;           // initial centre
    float const ycen0 = ycen;           //                of object

    double sigma11W = 1.5;              // quadratic moments of the
    double sigma12W = 0.0;              //     weighting fcn;
    double sigma22W = 1.5;              //               xx, xy, and yy

    double w11 = -1, w12 = -1, w22 = -1;        // current weights for moments; always set when iter == 0
    float e1_old = 1e6, e2_old = 1e6;           // old values of shape parameters e1 and e2
    float sigma11_ow_old = 1e6;                 // previous version of sigma11_ow

    typename ImageAdaptor<ImageT>::Image const &image = ImageAdaptor<ImageT>().getImage(mimage);

    if (std::isnan(xcen) || std::isnan(ycen)) {
        // Can't do anything
        shape->flags[SdssShapeAlgorithm::UNWEIGHTED_BAD.number] = true;
        return false;
    }

    bool interpflag = false;            // interpolate finer than a pixel?
    lsst::afw::geom::BoxI bbox;
    int iter = 0;                       // iteration number
    for (; iter < maxIter; iter++) {
        bbox = computeAdaptiveMomentsBBox(image.getBBox(afw::image::LOCAL), afw::geom::Point2D(xcen, ycen),
                                          sigma11W, sigma12W, sigma22W);
        std::tuple<std::pair<bool, double>, double, double, double> weights =
            getWeights(sigma11W, sigma12W, sigma22W);
        if (!std::get<0>(weights).first) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }

        double const detW = std::get<0>(weights).second;

#if 0                                   // this form was numerically unstable on my G4 powerbook
        assert(detW >= 0.0);
#else
        assert(sigma11W*sigma22W >= sigma12W*sigma12W - std::numeric_limits<float>::epsilon());
#endif

        {
            const double ow11 = w11;    // old
            const double ow12 = w12;    //     values
            const double ow22 = w22;    //            of w11, w12, w22

            w11 = std::get<1>(weights);
            w12 = std::get<2>(weights);
            w22 = std::get<3>(weights);

            if (shouldInterp(sigma11W, sigma22W, detW)) {
                if (!interpflag) {
                    interpflag = true;       // N.b.: stays set for this object
                    if (iter > 0) {
                        sigma11_ow_old = 1.e6; // force at least one more iteration
                        w11 = ow11;
                        w12 = ow12;
                        w22 = ow22;
                        iter--;         // we didn't update wXX
                    }
                }
            }
        }

        if (calcmom<false>(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22,
                           &I0, &sum, &sumx, &sumy, &sumxx, &sumxy, &sumyy, &sums4, negative) < 0) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }

#if 0
/*
 * Find new centre
 *
 * This is only needed if we update the centre; if we use the input position we've already done the work
 */
        xcen = sumx/sum;
        ycen = sumy/sum;
#endif
        shape->x = sumx/sum; // update centroid.  N.b. we're not setting errors here
        shape->y = sumy/sum;

        if (fabs(shape->x - xcen0) > shiftmax || fabs(shape->y - ycen0) > shiftmax) {
            shape->flags[SdssShapeAlgorithm::SHIFT.number] = true;
        }
/*
 * OK, we have the centre. Proceed to find the second moments.
 */
        float const sigma11_ow = sumxx/sum; // quadratic moments of
        float const sigma22_ow = sumyy/sum; //          weight*object
        float const sigma12_ow = sumxy/sum; //                 xx, xy, and yy

        if (sigma11_ow <= 0 || sigma22_ow <= 0) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }

        float const d = sigma11_ow + sigma22_ow; // current values of shape parameters
        float const e1 = (sigma11_ow - sigma22_ow)/d;
        float const e2 = 2.0*sigma12_ow/d;
/*
 * Did we converge?
 */
        if (iter > 0 &&
           fabs(e1 - e1_old) < tol1 && fabs(e2 - e2_old) < tol1 &&
           fabs(sigma11_ow/sigma11_ow_old - 1.0) < tol2 ) {
            break;                              // yes; we converged
        }

        e1_old = e1;
        e2_old = e2;
        sigma11_ow_old = sigma11_ow;
/*
 * Didn't converge, calculate new values for weighting function
 *
 * The product of two Gaussians is a Gaussian:
 * <x^2 exp(-a x^2 - 2bxy - cy^2) exp(-Ax^2 - 2Bxy - Cy^2)> =
 *                            <x^2 exp(-(a + A) x^2 - 2(b + B)xy - (c + C)y^2)>
 * i.e. the inverses of the covariances matrices add.
 *
 * We know sigmaXX_ow and sigmaXXW, the covariances of the weighted object
 * and of the weights themselves.  We can estimate the object's covariance as
 *   sigmaXX_ow^-1 - sigmaXXW^-1
 * and, as we want to find a set of weights with the _same_ covariance as the
 * object we take this to be the an estimate of our correct weights.
 *
 * N.b. This assumes that the object is roughly Gaussian.
 * Consider the object:
 *   O == delta(x + p) + delta(x - p)
 * the covariance of the weighted object is equal to that of the unweighted
 * object, and this prescription fails badly.  If we detect this, we set
 * the UNWEIGHTED flag, and calculate the UNweighted moments
 * instead.
 */
        {
            float n11, n12, n22;                // elements of inverse of next guess at weighting function
            float ow11, ow12, ow22;             // elements of inverse of sigmaXX_ow

            std::tuple<std::pair<bool, double>, double, double, double> weights =
                getWeights(sigma11_ow, sigma12_ow, sigma22_ow);
            if (!std::get<0>(weights).first) {
                shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
                break;
            }

            ow11 = std::get<1>(weights);
            ow12 = std::get<2>(weights);
            ow22 = std::get<3>(weights);

            n11 = ow11 - w11;
            n12 = ow12 - w12;
            n22 = ow22 - w22;

            weights = getWeights(n11, n12, n22);
            if (!std::get<0>(weights).first) {
                // product-of-Gaussians assumption failed
                shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
                break;
            }

            sigma11W = std::get<1>(weights);
            sigma12W = std::get<2>(weights);
            sigma22W = std::get<3>(weights);
        }

        if (sigma11W <= 0 || sigma22W <= 0) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
            break;
        }
    }

    if (iter == maxIter) {
        shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
        shape->flags[SdssShapeAlgorithm::MAXITER.number] = true;
    }

    if (sumxx + sumyy == 0.0) {
        shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = true;
    }
/*
 * Problems; try calculating the un-weighted moments
 */
    if (shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number]) {
        w11 = w22 = w12 = 0;
        if (calcmom<false>(image, xcen, ycen, bbox, bkgd, interpflag, w11, w12, w22,
                           &I0, &sum, &sumx, &sumy, &sumxx, &sumxy, &sumyy, NULL, negative) < 0 ||
	    (!negative && sum <= 0) || (negative && sum >= 0)) {
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number] = false;
            shape->flags[SdssShapeAlgorithm::UNWEIGHTED_BAD.number] = true;

            if (sum > 0) {
                shape->xx = 1/12.0;      // a single pixel
                shape->xy = 0.0;
                shape->yy = 1/12.0;
            }

            return false;
        }

        sigma11W = sumxx/sum;          // estimate of object moments
        sigma12W = sumxy/sum;          //   usually, object == weight
        sigma22W = sumyy/sum;          //      at this point
    }

    shape->flux = I0;
    shape->xx = sigma11W;
    shape->xy = sigma12W;
    shape->yy = sigma22W;

    if (shape->xx + shape->yy != 0.0) {
        int const ix = lsst::afw::image::positionToIndex(xcen);
        int const iy = lsst::afw::image::positionToIndex(ycen);

        if (ix >= 0 && ix < mimage.getWidth() && iy >= 0 && iy < mimage.getHeight()) {
            float const bkgd_var =
                ImageAdaptor<ImageT>().getVariance(mimage, ix, iy); // XXX Overestimate as it includes object

            if (bkgd_var > 0.0) {                                   // NaN is not > 0.0
                if (!(shape->flags[SdssShapeAlgorithm::UNWEIGHTED.number])) {
                    Matrix4d fisher = calc_fisher(*shape, bkgd_var); // Fisher matrix
                    Matrix4d cov = fisher.inverse();
                    // convention in afw::geom::ellipses is to order moments (xx, yy, xy),
                    // but the older algorithmic code uses (xx, xy, yy) - the order of
                    // indices here is not a bug.
                    shape->fluxSigma = std::sqrt(cov(0, 0));
                    shape->xxSigma = std::sqrt(cov(1, 1));
                    shape->xySigma = std::sqrt(cov(2, 2));
                    shape->yySigma = std::sqrt(cov(3, 3));
                    shape->flux_xx_Cov = cov(0, 1);
                    shape->flux_xy_Cov = cov(0, 2);
                    shape->flux_yy_Cov = cov(0, 3);
                    shape->xx_yy_Cov = cov(1, 3);
                    shape->xx_xy_Cov = cov(1, 2);
                    shape->yy_xy_Cov = cov(2, 3);
                }
            }
        }
    }

    return true;
}

} // anonymous


SdssShapeResult::SdssShapeResult() :
    flux_xx_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    flux_yy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    flux_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}


SdssShapeResultKey SdssShapeResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    bool doMeasurePsf
) {

    SdssShapeResultKey r;
    r._shapeResult = ShapeResultKey::addFields(schema, name, "elliptical Gaussian adaptive moments",
                                               SIGMA_ONLY);
    r._centroidResult = CentroidResultKey::addFields(schema, name, "elliptical Gaussian adaptive moments",
                                                     NO_UNCERTAINTY);
    r._fluxResult = FluxResultKey::addFields(schema, name, "elliptical Gaussian adaptive moments");

    // Only include PSF shape fields if doMeasurePsf = True
    if (doMeasurePsf) {
        r._includePsf = true;
        r._psfShapeResult = afwTable::QuadrupoleKey::addFields(
            schema, schema.join(name, "psf"), "adaptive moments of the PSF model at the object position");
    } else {
        r._includePsf = false;
    }

    r._flux_xx_Cov = schema.addField<ErrElement>(
        schema.join(name, "flux", "xx", "Cov"),
        (boost::format("uncertainty covariance between %s and %s")
         % schema.join(name, "flux") % schema.join(name, "xx")).str(),
        "count*pixel^2"
    );
    r._flux_yy_Cov = schema.addField<ErrElement>(
        schema.join(name, "flux", "yy", "Cov"),
        (boost::format("uncertainty covariance between %s and %s")
         % schema.join(name, "flux") % schema.join(name, "yy")).str(),
        "count*pixel^2"
    );
    r._flux_xy_Cov = schema.addField<ErrElement>(
        schema.join(name, "flux", "xy", "Cov"),
        (boost::format("uncertainty covariance between %s and %s")
         % schema.join(name, "flux") % schema.join(name, "xy")).str(),
        "count*pixel^2"
    );

    // Skip the psf flag if not recording the PSF shape.
    if (r._includePsf) {
        r._flagHandler = FlagHandler::addFields(schema, name, SdssShapeAlgorithm::getFlagDefinitions());
    }
    else {
        r._flagHandler = FlagHandler::addFields(schema, name, SdssShapeAlgorithm::getFlagDefinitions(),
                            {SdssShapeAlgorithm::PSF_SHAPE_BAD});
    }
    return r;
}

SdssShapeResultKey::SdssShapeResultKey(afw::table::SubSchema const & s) :
    _shapeResult(s),
    _centroidResult(s),
    _fluxResult(s),
    _flux_xx_Cov(s["flux"]["xx"]["Cov"]),
    _flux_yy_Cov(s["flux"]["yy"]["Cov"]),
    _flux_xy_Cov(s["flux"]["xy"]["Cov"])
{
    // The input SubSchema may optionally provide for a PSF.
    try {
        _psfShapeResult = afwTable::QuadrupoleKey(s["psf"]);
        _flagHandler = FlagHandler(s, SdssShapeAlgorithm::getFlagDefinitions());
        _includePsf = true;
    } catch (pex::exceptions::NotFoundError& e) {
        _flagHandler = FlagHandler(s, SdssShapeAlgorithm::getFlagDefinitions(), {SdssShapeAlgorithm::PSF_SHAPE_BAD});
        _includePsf = false;
    }
}

SdssShapeResult SdssShapeResultKey::get(afw::table::BaseRecord const & record) const {
    SdssShapeResult result;
    static_cast<ShapeResult&>(result) = record.get(_shapeResult);
    static_cast<CentroidResult&>(result) = record.get(_centroidResult);
    static_cast<FluxResult&>(result) = record.get(_fluxResult);
    result.flux_xx_Cov = record.get(_flux_xx_Cov);
    result.flux_yy_Cov = record.get(_flux_yy_Cov);
    result.flux_xy_Cov = record.get(_flux_xy_Cov);
    for (unsigned int n = 0; n < SdssShapeAlgorithm::N_FLAGS; ++n) {
        if (n == SdssShapeAlgorithm::PSF_SHAPE_BAD.number && !_includePsf) continue;
        result.flags[n] = _flagHandler.getValue(record, n);
    }
    return result;
}

afwGeom::ellipses::Quadrupole SdssShapeResultKey::getPsfShape(afw::table::BaseRecord const & record) const {
    return record.get(_psfShapeResult);
}

void SdssShapeResultKey::set(afw::table::BaseRecord & record, SdssShapeResult const & value) const {
    record.set(_shapeResult, value);
    record.set(_centroidResult, value);
    record.set(_fluxResult, value);
    record.set(_flux_xx_Cov, value.flux_xx_Cov);
    record.set(_flux_yy_Cov, value.flux_yy_Cov);
    record.set(_flux_xy_Cov, value.flux_xy_Cov);
    for (unsigned int n = 0; n < SdssShapeAlgorithm::N_FLAGS; ++n) {
        if (n == SdssShapeAlgorithm::PSF_SHAPE_BAD.number && !_includePsf) continue;
        _flagHandler.setValue(record, n, value.flags[n]);
    }
}

void SdssShapeResultKey::setPsfShape(afw::table::BaseRecord & record,
                                     afwGeom::ellipses::Quadrupole const & value) const {
    record.set(_psfShapeResult, value);
}

bool SdssShapeResultKey::operator==(SdssShapeResultKey const & other) const {
    return _shapeResult == other._shapeResult &&
        _centroidResult == other._centroidResult &&
        _fluxResult == other._fluxResult &&
        _psfShapeResult == other._psfShapeResult &&
        _flux_xx_Cov == other._flux_xx_Cov &&
        _flux_yy_Cov == other._flux_yy_Cov &&
        _flux_xy_Cov == other._flux_xy_Cov;
    // don't bother with flags - if we've gotten this far, it's basically impossible the flags don't match
}

bool SdssShapeResultKey::isValid() const {
    return _shapeResult.isValid() &&
        _centroidResult.isValid() &&
        _fluxResult.isValid() &&
        _psfShapeResult.isValid() &&
        _flux_xx_Cov.isValid() &&
        _flux_yy_Cov.isValid() &&
        _flux_xy_Cov.isValid();
    // don't bother with flags - if we've gotten this far, it's basically impossible the flags are invalid
}

SdssShapeAlgorithm::SdssShapeAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
)
  : _ctrl(ctrl),
    _resultKey(ResultKey::addFields(schema, name, ctrl.doMeasurePsf)),
    _centroidExtractor(schema, name)
{}

template <typename ImageT>
SdssShapeResult SdssShapeAlgorithm::computeAdaptiveMoments(
    ImageT const & image,
    afw::geom::Point2D const & center,
    bool negative,
    Control const & control
) {
    double xcen = center.getX();         // object's column position
    double ycen = center.getY();         // object's row position

    xcen -= image.getX0();               // work in image Pixel coordinates
    ycen -= image.getY0();

    float shiftmax = control.maxShift;   // Max allowed centroid shift
    if (shiftmax < 2) {
        shiftmax = 2;
    } else if (shiftmax > 10) {
        shiftmax = 10;
    }

    SdssShapeResult result;
    try {
        result.flags[FAILURE.number] = !getAdaptiveMoments(
            image, control.background, xcen, ycen, shiftmax, &result,
            control.maxIter, control.tol1, control.tol2, negative
        );
    } catch (pex::exceptions::Exception & err) {
        result.flags[FAILURE.number] = true;
    }
    if (result.flags[UNWEIGHTED.number] || result.flags[SHIFT.number]) {
        // These are also considered fatal errors in terms of the quality of the results,
        // even though they do produce some results.
        result.flags[FAILURE.number] = true;
    }
    if (result.getQuadrupole().getIxx()*result.getQuadrupole().getIyy() <
            (1.0 + 1.0e-6)*result.getQuadrupole().getIxy()*result.getQuadrupole().getIxy())
                  // We are checking that Ixx*Iyy > (1 + epsilon)*Ixy*Ixy where epsilon is suitably small. The
                  // value of epsilon used here is a magic number. DM-5801 is supposed to figure out if we are
                  // to keep this value.
        {
        if (!result.flags[FAILURE.number]) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicError,
                "Should not get singular moments unless a flag is set");
        }
    }

    // getAdaptiveMoments() just computes the zeroth moment in result.flux (and its error in
    // result.fluxSigma, result.flux_xx_Cov, etc.)  That's related to the flux by some geometric
    // factors, which we apply here.
    double fluxScale = computeFluxScale(result);

    result.flux *= fluxScale;
    result.fluxSigma *= fluxScale;
    result.x += image.getX0();
    result.y += image.getY0();

    if (ImageAdaptor<ImageT>::hasVariance) {
        result.flux_xx_Cov *= fluxScale;
        result.flux_yy_Cov *= fluxScale;
        result.flux_xy_Cov *= fluxScale;
    }

    return result;
}

template <typename ImageT>
FluxResult SdssShapeAlgorithm::computeFixedMomentsFlux(
    ImageT const & image,
    afw::geom::ellipses::Quadrupole const & shape,
    afw::geom::Point2D const & center
) {
    // while arguments to computeFixedMomentsFlux are in PARENT coordinates, the implementation is LOCAL.
    afw::geom::Point2D localCenter = center - afw::geom::Extent2D(image.getXY0());

    afwGeom::BoxI const bbox = computeAdaptiveMomentsBBox(image.getBBox(afw::image::LOCAL),
                                                          localCenter,
                                                          shape.getIxx(), shape.getIxy(), shape.getIyy());

    std::tuple<std::pair<bool, double>, double, double, double> weights =
        getWeights(shape.getIxx(), shape.getIxy(), shape.getIyy());

    FluxResult result;

    if (!std::get<0>(weights).first) {
        throw pex::exceptions::InvalidParameterError("Input shape is singular");
    }

    double const w11 = std::get<1>(weights);
    double const w12 = std::get<2>(weights);
    double const w22 = std::get<3>(weights);
    bool const interp = shouldInterp(shape.getIxx(), shape.getIyy(), std::get<0>(weights).second);

    double i0 = 0;                      // amplitude of Gaussian
    if (calcmom<true>(ImageAdaptor<ImageT>().getImage(image), localCenter.getX(), localCenter.getY(),
                      bbox, 0.0, interp, w11, w12, w22, &i0, NULL, NULL, NULL, NULL, NULL, NULL, NULL)< 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, "Error from calcmom");
    }

    double const wArea = afw::geom::PI*std::sqrt(shape.getDeterminant());

    result.flux = i0*2*wArea;

    if (ImageAdaptor<ImageT>::hasVariance) {
        int ix = static_cast<int>(center.getX() - image.getX0());
        int iy = static_cast<int>(center.getY() - image.getY0());
        if (!image.getBBox(afw::image::LOCAL).contains(afw::geom::Point2I(ix, iy))) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                              (boost::format("Center (%d,%d) not in image (%dx%d)") %
                               ix % iy % image.getWidth() % image.getHeight()).str());
        }
        double var = ImageAdaptor<ImageT>().getVariance(image, ix, iy);
        double i0Err = std::sqrt(var/wArea);
        result.fluxSigma = i0Err*2*wArea;
    }

    return result;
}

void SdssShapeAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    bool negative = false;

    try {
        negative = measRecord.get(measRecord.getSchema().find<afw::table::Flag>("flags_negative").key);
    } catch(pexExcept::Exception &e) {
    }
    SdssShapeResult result = computeAdaptiveMoments(
        exposure.getMaskedImage(),
        _centroidExtractor(measRecord, _resultKey.getFlagHandler()),
        negative,
        _ctrl
    );

    if (_ctrl.doMeasurePsf) {
        // Compute moments of Psf model.  In the interest of implementing this quickly, we're just
        // calling Psf::computeShape(), which delegates to SdssShapeResult::computeAdaptiveMoments
        // for all nontrivial Psf classes.  But this could in theory save the results of a shape
        // computed some other way as part of base_SdssShape, which might be confusing.  We should
        // fix this eventually either by making Psf shape measurement not part of base_SdssShape, or
        // by making the measurements stored with shape.sdss always computed via the
        // SdssShapeAlgorithm instead of delegating to the Psf class.
        try {
            PTR(afw::detection::Psf const) psf = exposure.getPsf();
            if (!psf) {
                result.flags[PSF_SHAPE_BAD.number] = true;
            } else {
                _resultKey.setPsfShape(measRecord, psf->computeShape(afw::geom::Point2D(result.x, result.y)));
            }
        } catch (pex::exceptions::Exception & err) {
            result.flags[PSF_SHAPE_BAD.number] = true;
        }
    }

    measRecord.set(_resultKey, result);
}

void SdssShapeAlgorithm::fail(
    afw::table::SourceRecord & measRecord,
    MeasurementError * error
) const {
    _resultKey.getFlagHandler().handleFailure(measRecord, error);
}

#define INSTANTIATE_IMAGE(IMAGE) \
    template SdssShapeResult SdssShapeAlgorithm::computeAdaptiveMoments( \
        IMAGE const &,                                                  \
        afw::geom::Point2D const &,                                     \
        bool,                                                           \
        Control const &                                                 \
    );                                                                  \
    template FluxResult SdssShapeAlgorithm::computeFixedMomentsFlux( \
        IMAGE const &,                                                  \
        afw::geom::ellipses::Quadrupole const &,                        \
        afw::geom::Point2D const &                                      \
    )

#define INSTANTIATE_PIXEL(PIXEL) \
    INSTANTIATE_IMAGE(lsst::afw::image::Image<PIXEL>); \
    INSTANTIATE_IMAGE(lsst::afw::image::MaskedImage<PIXEL>);

INSTANTIATE_PIXEL(int);
INSTANTIATE_PIXEL(float);
INSTANTIATE_PIXEL(double);

SdssShapeTransform::SdssShapeTransform(
    Control const & ctrl,
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    BaseTransform{name},
    _fluxTransform{name, mapper},
    _centroidTransform{name, mapper}
{
    // If the input schema has a PSF flag -- it's optional --  assume we are also transforming the PSF.
    _transformPsf = mapper.getInputSchema().getNames().count("sdssShape_flag_psf") ? true : false;

    // Skip the last flag if not transforming the PSF shape.
    for (std::size_t i = 0; i < SdssShapeAlgorithm::getFlagDefinitions().size(); i++) {
        FlagDefinition const & flag = SdssShapeAlgorithm::getFlagDefinitions()[i];
        if (flag == SdssShapeAlgorithm::FAILURE) continue;
        if (mapper.getInputSchema().getNames().count(name + "_" + flag.name) == 0) continue;
        afw::table::Key<afw::table::Flag> key = mapper.getInputSchema().find<afw::table::Flag>(
            name + "_" + flag.name).key;
        mapper.addMapping(key);
    }

    _outShapeKey = ShapeResultKey::addFields(mapper.editOutputSchema(), name, "Shape in celestial moments",
                                             FULL_COVARIANCE, afw::table::CoordinateType::CELESTIAL);
    if (_transformPsf) {
        _outPsfShapeKey = afwTable::QuadrupoleKey::addFields(mapper.editOutputSchema(), name + "_psf",
                                                             "PSF shape in celestial moments",
                                                             afw::table::CoordinateType::CELESTIAL);
    }
}

void SdssShapeTransform::operator()(
    afw::table::SourceCatalog const & inputCatalog,
    afw::table::BaseCatalog & outputCatalog,
    afw::image::Wcs const & wcs,
    afw::image::Calib const & calib
) const {
    // The flux and cetroid transforms will check that the catalog lengths
    // match and throw if not, so we don't repeat the test here.
    _fluxTransform(inputCatalog, outputCatalog, wcs, calib);
    _centroidTransform(inputCatalog, outputCatalog, wcs, calib);

    CentroidResultKey centroidKey(inputCatalog.getSchema()[_name]);
    ShapeResultKey inShapeKey(inputCatalog.getSchema()[_name]);
    afwTable::QuadrupoleKey inPsfShapeKey;
    if (_transformPsf) {
        inPsfShapeKey = afwTable::QuadrupoleKey(
            inputCatalog.getSchema()[inputCatalog.getSchema().join(_name, "psf")]);
    }

    afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
    afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();
    for (; inSrc != inputCatalog.end(); ++inSrc, ++outSrc) {
        ShapeResult inShape = inShapeKey.get(*inSrc);
        ShapeResult outShape;

        // The transformation from the (x, y) to the (Ra, Dec) basis.
        afw::geom::AffineTransform crdTr = wcs.linearizePixelToSky(centroidKey.get(*inSrc).getCentroid(),
                                                                   afw::geom::radians);
        outShape.setShape(inShape.getShape().transform(crdTr.getLinear()));

        // Transformation matrix from pixel to celestial basis.
        ShapeTrMatrix m = makeShapeTransformMatrix(crdTr.getLinear());
        outShape.setShapeErr((m*inShape.getShapeErr().cast<double>()*m.transpose()).cast<ErrElement>());

        _outShapeKey.set(*outSrc, outShape);

        if (_transformPsf) {
            _outPsfShapeKey.set(*outSrc, inPsfShapeKey.get(*inSrc).transform(crdTr.getLinear()));
        }
    }
}

}}} // end namespace lsst::meas::base
