// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST.
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

#ifndef LSST_MEAS_BASE_constants_h_INCLUDED
#define LSST_MEAS_BASE_constants_h_INCLUDED

#include "Eigen/Core"

#include "lsst/pex/exceptions.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief An enum used to specify how much uncertainty information measurement algorithms provide.
 *
 *  Currently, only ResultMappers (not Results) make use of these distinctions; Result structs always
 *  have data members that could hold the full-covariance, but may set some of these to NaN.
 */
enum UncertaintyEnum {
    NO_UNCERTAINTY = 0, ///< Algorithm provides no uncertainy information at all
    SIGMA_ONLY = 1,     ///< Only the diagonal elements of the covariance matrix are provided
    FULL_COVARIANCE = 2 ///< The full covariance matrix is provided
};

//@{ Typedefs that define the C++ types we typically use for common measurements
typedef int ElementCount;
typedef double Flux;
typedef double FluxErrElement;
typedef double Mag;
typedef double MagErrElement;
typedef float ErrElement;
typedef double CentroidElement;
typedef double ShapeElement;
typedef geom::Point<CentroidElement,2> Centroid;
typedef Eigen::Matrix<ErrElement,2,2,Eigen::DontAlign> CentroidCov;
typedef afw::geom::ellipses::Quadrupole Shape;
typedef Eigen::Matrix<ErrElement,3,3,Eigen::DontAlign> ShapeCov;
typedef Eigen::Matrix<ShapeElement,3,3,Eigen::DontAlign> ShapeTrMatrix;
//@}

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_constants_h_INCLUDED
