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
#include "lsst/meas/base/exceptions.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief Empty control object used by algorithm classes that don't have any configuration parameters.
 *
 *  It would be a bit cleaner in C++ just to not have one, but having a null one makes the Python side
 *  much cleaner.
 */
struct NullControl {};

//@{
/**
 *  These reusable structs represent predefined Inputs which will work for most algorithms.
 *
 *  Each Input can be constructed directly from its constituents, or from a SourceRecord.
 *  The latter lets us use these interchangeably in Python.  They all also provide a static
 *  method to create a std::vector of Inputs from a SourceCatalog.
 */

/// An Input struct for algorithms that require only a Footprint
struct FootprintInput {
    typedef std::vector<FootprintInput> Vector;

    PTR(afw::detection::Footprint) footprint;

    explicit FootprintInput(PTR(afw::detection::Footprint) footprint_) : footprint(footprint_) {}

    explicit FootprintInput(afw::table::SourceRecord const & record) : footprint(record.getFootprint()) {}

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

    bool hasFlaggedDependencies() const { return false; }

};

/// An Input struct for algorithms that require a position as well as a Footprint
struct FootprintCentroidInput : public FootprintInput {
    typedef std::vector<FootprintCentroidInput> Vector;

    afw::geom::Point2D position;

    explicit FootprintCentroidInput(
        PTR(afw::detection::Footprint) footprint_, afw::geom::Point2D const & position_
    ) : FootprintInput(footprint_), position(position_) {}

    FootprintCentroidInput(afw::table::SourceRecord const & record);

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

    bool hasFlaggedDependencies() const { return false; }

};

/// An Input struct for algorithms that require a position and shape as well as a Footprint
struct FootprintCentroidShapeInput : public FootprintCentroidInput {
    typedef std::vector<FootprintCentroidShapeInput> Vector;

    afw::geom::ellipses::Quadrupole shape;
    bool shapeFlag;

    FootprintCentroidShapeInput(
        PTR(afw::detection::Footprint) footprint_,
        afw::geom::Point2D const & position_,
        afw::geom::ellipses::Quadrupole const & shape_,
        bool const & shapeFlag_
    ) : FootprintCentroidInput(footprint_, position_), shape(shape_), shapeFlag(shapeFlag_) {}

    explicit FootprintCentroidShapeInput(afw::table::SourceRecord const & record) :
        FootprintCentroidInput(record)
    {
        if (record.getTable()->hasShapeSlot()) {
            shape = record.getShape();
            shapeFlag = record.getShapeFlag();
        } else {
            throw LSST_EXCEPT(
                lsst::meas::base::FatalAlgorithmError,
                "Shape slot must be setup for algorithm to be run"
            );
        }
    }

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

    bool hasFlaggedDependencies() const { return shapeFlag; }

};

//@}

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Inputs_h_INCLUDED
