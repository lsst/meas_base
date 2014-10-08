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

#ifndef LSST_MEAS_BASE_ShapeUtilities_h_INCLUDED
#define LSST_MEAS_BASE_ShapeUtilities_h_INCLUDED

#include "lsst/meas/base/constants.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A reusable component for result structs for moments-based shape measurements.
 *
 *  Shape measurements and their errors should always be in pixels coordinates.  This struct should generally
 *  be preferred over a custom struct with other ellipse parametrizations unless the measurement takes place
 *  in another parametrization and a transformation to this one would result in a loss of information or
 *  obfuscate the results of the measurement (i.e. use this one unless you have a good reason not to).
 */
struct ShapeResult {
    ShapeElement xx; // image or model second moment for x^2
    ShapeElement yy; // image or model second moment for y^2
    ShapeElement xy; // image or model second moment for xy^2
    ErrElement xxSigma; ///< 1-Sigma uncertainty on xx (sqrt of variance)
    ErrElement yySigma; ///< 1-Sigma uncertainty on yy (sqrt of variance)
    ErrElement xySigma; ///< 1-Sigma uncertainty on xy (sqrt of variance)
    ErrElement xx_yy_Cov; ///< xx,yy term in the uncertainty convariance matrix
    ErrElement xx_xy_Cov; ///< xx,xy term in the uncertainty convariance matrix
    ErrElement yy_xy_Cov; ///< yy,xy term in the uncertainty convariance matrix

    /**
     *  @brief Return an afw::geom::ellipses object corresponding to xx, yy, xy.
     *
     *  This method can be used to return an average radius for the measured shape, e.g.
     *  @c getShape().getDeterminantRadius()
     */
    Shape const getShape() const { return Shape(xx, yy, xy); }

    /// Return the 3x3 symmetric covariance matrix, with rows and columns ordered (xx, yy, xy)
    ShapeCov const getCov() const {
        ShapeCov m;
        m <<
            xxSigma*xxSigma, xx_yy_Cov, xx_xy_Cov,
            xx_yy_Cov, yySigma*yySigma, yy_xy_Cov,
            xx_xy_Cov, yy_xy_Cov, xySigma*xySigma;
        return m;
    }

    ShapeResult(); ///< Constructor; initializes everything to NaN.

};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ShapeUtilities_h_INCLUDED
