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

#include <bitset>

#include "boost/array.hpp"
#include "Eigen/Core"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"

namespace lsst { namespace meas { namespace base {

class MeasurementError : public pex::exceptions::RuntimeErrorException {
public:
    MeasurementError(LSST_EARGS_TYPED, std::size_t flagBit) :
        pex::exceptions::RuntimeErrorException(LSST_EARGS_UNTYPED),
        _flagBit(flagBit)
    {}

    std::size_t getFlagBit() const { return _flagBit; }

    virtual char const* getType(void) const throw() { return "lsst::meas::base::MeasurementError *"; };

    virtual lsst::pex::exceptions::Exception* clone(void) const {
        return new MeasurementError(*this);
    };
    
private:
    std::size_t _flagBit;
};


// Typedefs that define the C++ types we typically use for common measurements

typedef float Flux;
typedef float ErrElement;
typedef double CentroidElement;
typedef double ShapeElement;
typedef afw::geom::Point<CentroidElement,2> Centroid;
typedef Eigen::Matrix<ErrElement,2,2,Eigen::DontAlign> CentroidCov;
typedef afw::geom::ellipses::Quadrupole Shape;
typedef Eigen::Matrix<ErrElement,3,3,Eigen::DontAlign> ShapeCov;

// We expect the structs below will be reused (used directly, subclassed, or composition) by most algorithms.
// In the measurement framework, each algorithm should also have at least one flag field, but that will be
// added by the plugin wrapper layer and set when the algorithm code here throws an exception.

// Simple POD struct used to define flags for use by the ResultMapper classes later.
// When we switch to C++11, we can make these std::string, but for now we'll use
// C-strings so we can initialize arrays of these using initializer lists.
struct FlagDef {
    char const * name;
    char const * doc;
};

template <std::size_t N>
struct FlagsResult {
    std::bitset<N> flags;
};

struct FluxResult {
    Flux flux;
    ErrElement fluxSigma;

    FluxResult();
};

struct CentroidResult {
    CentroidElement x;
    CentroidElement y;
    ErrElement xSigma;
    ErrElement ySigma;
    ErrElement x_y_Cov;

    Centroid const getCentroid() const { return Centroid(x, y); }

    CentroidCov const getCov() const {
        CentroidCov m;
        m <<
            xSigma*xSigma, x_y_Cov,
            x_y_Cov, ySigma*ySigma;
        return m;
    }

    CentroidResult();

};

struct ShapeResult {
    ShapeElement xx;
    ShapeElement yy;
    ShapeElement xy;
    ErrElement xxSigma;
    ErrElement yySigma;
    ErrElement xySigma;
    ErrElement xx_yy_Cov;
    ErrElement xx_xy_Cov;
    ErrElement yy_xy_Cov;

    Shape const getShape() const { return Shape(xx, yy, xy); }

    ShapeCov const getCov() const {
        ShapeCov m;
        m <<
            xxSigma*xxSigma, xx_yy_Cov, xx_xy_Cov,
            xx_yy_Cov, yySigma*yySigma, yy_xy_Cov,
            xx_xy_Cov, yy_xy_Cov, xySigma*xySigma;
        return m;
    }

    ShapeResult();

};

template <typename Algorithm, typename T1>
struct SimpleResult1 : public T1, public FlagsResult<Algorithm::N_FLAGS> {};

template <typename Algorithm, typename T1, typename T2>
struct SimpleResult2 : public T1, public T2, public FlagsResult<Algorithm::N_FLAGS> {};

template <typename Algorithm, typename T1, typename T2, typename T3>
struct SimpleResult3 : public T1, public T2, public T3, public FlagsResult<Algorithm::N_FLAGS> {};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Results_h_INCLUDED
