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

#ifndef LSST_MEAS_BASE_Inputs_h_INCLUDED
#define LSST_MEAS_BASE_Inputs_h_INCLUDED

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst { namespace meas { namespace base {

// Reusable structs that represent the inputs for many algorithms, used mainly so we can build
// a vector of inputs for multi-object measurement cases.

struct AlgorithmInput1 {
    PTR(afw::detection::Footprint) footprint;
};

struct AlgorithmInput2 : public AlgorithmInput1 {
    PTR(afw::detection::Footprint) footprint;
    afw::geom::Point2D position;
};

struct AlgorithmInput3 : public AlgorithmInput2 {
    afw::geom::ellipses::Quadrupole shape;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Inputs_h_INCLUDED
