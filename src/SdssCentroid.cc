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
#include <iostream>
#include <cmath>
#include <numeric>
#include "ndarray/eigen.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/math/ConvolveImage.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/base/SdssCentroid.h"
#include "lsst/meas/base/algorithms/SdssCentroidTemplates.h"


namespace lsst { namespace meas { namespace base {


SdssCentroidAlgorithm::ResultMapper SdssCentroidAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, FULL_COVARIANCE);
}

template <typename T>
SdssCentroidAlgorithm::Result SdssCentroidAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Control const & ctrl
) {
    Result result;

    // This code has been moved essentially without change from meas_algorithms
    // The only changes were:
    // Change the exceptions to MeasurementErrors with the correct flag bits
    // Change to set values in result rather than in a source record.
    result.x = center.getX(); result.y = center.getY(); // better than NaN
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    typedef typename MaskedImageT::Variance VarianceT;

    MaskedImageT const& mimage = exposure.getMaskedImage();
    ImageT const& image = *mimage.getImage();
    CONST_PTR(lsst::afw::detection::Psf) psf = exposure.getPsf();

    int const x = static_cast<int>(center.getX() + 0.5) - image.getX0(); // in image Pixel coords
    int const y = static_cast<int>(center.getY() + 0.5) - image.getY0();

    if (x < 0 || x >= image.getWidth() || y < 0 || y >= image.getHeight()) {
            throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[EDGE].doc,
            EDGE
        );
    }
    /*
     * If a PSF is provided, smooth the object with that PSF
     */
    if (!psf) {                  // image is presumably already smoothed
        // FIXME: the above logic is probably bad; this option should probably be a config parameter
        //psf.reset(new meas_algorithms::DoubleGaussianPsf(11, 11, 0.01));
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_PSF].doc,
            NO_PSF
        );
    }
    
    int binX = 1;
    int binY = 1;
    double xc=0., yc=0., dxc=0., dyc=0.;            // estimated centre and error therein
    for(int binsize = 1; binsize <= ctrl.binmax; binsize *= 2) {
        std::pair<MaskedImageT, double> result = algorithms::smoothAndBinImage(psf, x, y, mimage, binX, binY);
        MaskedImageT const smoothedImage = result.first;
        double const smoothingSigma = result.second;

        typename MaskedImageT::xy_locator mim = smoothedImage.xy_at(smoothedImage.getWidth()/2,
                                                                    smoothedImage.getHeight()/2);

        try {
            double sizeX2, sizeY2;      // object widths^2 in x and y directions
            double peakVal;             // peak intensity in image

            algorithms::doMeasureCentroidImpl(&xc, &dxc, &yc, &dyc, &sizeX2, &sizeY2, &peakVal, mim, smoothingSigma);

            if(binsize > 1) {
                // dilate from the lower left corner of central pixel
                xc = (xc + 0.5)*binX - 0.5;
                dxc *= binX;
                sizeX2 *= binX*binX;

                yc = (yc + 0.5)*binY - 0.5;
                dyc *= binY;
                sizeY2 *= binY*binY;
            }
      
            xc += x;                    // xc, yc are measured relative to pixel (x, y)
            yc += y;

            double const fac = ctrl.wfac*(1 + smoothingSigma*smoothingSigma);
            double const facX2 = fac*binX*binX;
            double const facY2 = fac*binY*binY;

            if (sizeX2 < facX2 && ::pow(xc - x, 2) < facX2 &&
                sizeY2 < facY2 && ::pow(yc - y, 2) < facY2) {
                if (binsize > 1 || ctrl.peakMin < 0.0 || peakVal > ctrl.peakMin) {
                    break;
                }
            }

            if (sizeX2 >= facX2 || ::pow(xc - x, 2) >= facX2) {
                binX *= 2;
            }
            if (sizeY2 >= facY2 || ::pow(yc - y, 2) >= facY2) {
                binY *= 2;
            }
        }
        catch(lsst::pex::exceptions::Exception &e) {
            throw LSST_EXCEPT(
                MeasurementError,
                getFlagDefinitions()[BAD_DATA].doc,
                BAD_DATA
            );
        }
    }
    result.x = lsst::afw::image::indexToPosition(xc + image.getX0());
    result.y = lsst::afw::image::indexToPosition(yc + image.getY0());

    result.xSigma = sqrt(dxc*dxc);
    result.ySigma = sqrt(dyc*dyc);
    // FIXME: should include off-diagonal term in covariance
    result.x_y_Cov = 0.0;
    return result;
}

template <typename T>
SdssCentroidAlgorithm::Result SdssCentroidAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template SdssCentroidAlgorithm::Result SdssCentroidAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    SdssCentroidAlgorithm::Result SdssCentroidAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

