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

#ifndef LSST_MEAS_BASE_Results_h_INCLUDED
#define LSST_MEAS_BASE_Results_h_INCLUDED

#include "Eigen/Core"

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"

namespace lsst { namespace meas { namespace base {

// Typedefs that define the C++ types we typically use for common measurements

typedef float Flux;
typedef float FluxErr;
typedef afw::geom::Point2D Centroid;
typedef Eigen::Matrix<float,2,2,Eigen::DontAlign> CentroidCov;
typedef afw::geom::ellipses::Quadrupole Shape;
typedef Eigen::Matrix<float,3,3,Eigen::DontAlign> ShapeCov;

// We expect the structs below will be reused (used directly, subclassed, or composition) by most algorithms.
// In the measurement framework, each algorithm should also have at least one flag field, but that will be
// added by the plugin wrapper layer and set when the algorithm code here throws an exception.

struct FluxAlgorithmResult {
    Flux value;
    FluxErr err;

    FluxAlgorithmResult(Flux value_, FluxErr err_) : value(value_), err(err_) {}

    FluxAlgorithmResult();
};

struct CentroidAlgorithmResult {
    Centroid value;
    CentroidCov cov;

    CentroidAlgorithmResult(
        Centroid const & value_,
        CentroidCov const & cov_
    ) : value(value_), cov(cov_) {}

    CentroidAlgorithmResult();

};

struct ShapeAlgorithmResult {
    Shape value;
    ShapeCov cov;

    ShapeAlgorithmResult(
        Shape const & value_,
        ShapeCov const & cov_
    ) : value(value_), cov(cov_) {}

    ShapeAlgorithmResult();

};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Results_h_INCLUDED
