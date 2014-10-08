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

#ifndef LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED
#define LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED

#include "lsst/meas/base/constants.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A reusable component for result structs for centroid or other position measurements.
 *
 *  Centroid measurements and their errors should always be in pixels, relative to the image's xy0.
 */
struct CentroidResult {
    CentroidElement x; ///< x (column) coordinate of the measured position
    CentroidElement y; ///< y (row) coordinate of the measured position
    ErrElement xSigma; ///< 1-Sigma uncertainty on x (sqrt of variance)
    ErrElement ySigma; ///< 1-Sigma uncertainty on y (sqrt of variance)
    ErrElement x_y_Cov; ///< x,y term in the uncertainty convariance matrix

    /// Return a Point object containing the measured x and y
    Centroid const getCentroid() const { return Centroid(x, y); }

    /// Return the 2x2 symmetric covariance matrix, with rows and columns ordered (x, y)
    CentroidCov const getCov() const {
        CentroidCov m;
        m <<
            xSigma*xSigma, x_y_Cov,
            x_y_Cov, ySigma*ySigma;
        return m;
    }

    CentroidResult(); ///< Constructor; initializes everything to NaN.

};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED
