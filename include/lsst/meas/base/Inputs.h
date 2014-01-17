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

#include <vector>

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/table/Source.h"

namespace lsst { namespace meas { namespace base {

struct NullControl {
    // Empty control object used by algorithm classes that don't have any configuration parameters
    // (it'd be cleaner in C++ just to not have one, but having a null one makes the Python side
    // *much* cleaner).
};

// Reusable structs that represent the inputs for many algorithms, used mainly so we can build
// a vector of inputs for multi-object measurement cases.

struct AlgorithmInput1 {
    typedef std::vector<AlgorithmInput1> Vector;

    PTR(afw::detection::Footprint) footprint;

    explicit AlgorithmInput1(afw::table::SourceRecord const & record) : footprint(record.getFootprint()) {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

struct AlgorithmInput2 : public AlgorithmInput1 {
    typedef std::vector<AlgorithmInput2> Vector;

    afw::geom::Point2D position;

    explicit AlgorithmInput2(afw::table::SourceRecord const & record) :
        AlgorithmInput1(record), position(record.getCentroid())
    {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

struct AlgorithmInput3 : public AlgorithmInput2 {
    typedef std::vector<AlgorithmInput3> Vector;

    afw::geom::ellipses::Quadrupole shape;

    explicit AlgorithmInput3(afw::table::SourceRecord const & record) :
        AlgorithmInput2(record), shape(record.getShape())
    {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Inputs_h_INCLUDED
