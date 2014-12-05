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

#ifndef LSST_MEAS_BASE_SillyCentroid_h_INCLUDED
#define LSST_MEAS_BASE_SillyCentroid_h_INCLUDED

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "lsst/pex/config.h"
#include "ndarray/eigen.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace test { namespace foo { namespace bar {
/**
 *  @brief A C++ control class to handle SillyCentroidAlgorithm's configuration
 */
class SillyCentroidControl {
public:
    LSST_CONTROL_FIELD(param, double, "Value to subtrace from the centroid position");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    SillyCentroidControl() : param(0.0) {}
};

/**
 *  @brief A class that calculates a centroid as a simple unweighted first moment
 *         of the 3x3 region around a pixel.
 *
 *   A fixed background (set via config) may optionally be subtracted.
 *   This algorithm does not currently report an error, but it probably should.
 */

class SillyCentroidAlgorithm : public lsst::meas::base::SimpleAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum {
        FAILURE=lsst::meas::base::FlagHandler::FAILURE,
        NO_COUNTS,
        EDGE,
        N_FLAGS
    };

    typedef SillyCentroidControl Control;

    SillyCentroidAlgorithm(
        Control const & ctrl,
        std::string const & name,
        lsst::afw::table::Schema & schema
    ) : _ctrl(ctrl),
        _centroidKey(
            lsst::meas::base::CentroidResultKey::addFields(schema, name, "centroid from Silly Centroid algorithm", lsst::meas::base::SIGMA_ONLY)
        ),
        _centroidExtractor(schema, name, true)
    {
        static boost::array<lsst::meas::base::FlagDefinition,N_FLAGS> const flagDefs = {{
            {"flag", "general failure flag, set if anything went wrong"},
            {"flag_noCounts", "Object to be centroided has no counts"},
            {"flag_edge", "Object too close to edge"}
        }};
        _flagHandler = lsst::meas::base::FlagHandler::addFields(schema, name, flagDefs.begin(), flagDefs.end());
    }
    

private:

    // These are private so they doesn't shadow the other overloads in base classes;
    // we can still call it via the public method on the base class.  We could have
    // used a using declaration instead, but Swig had trouble with that here.

    Control _ctrl;

    void measure(
        lsst::afw::table::SourceRecord & measRecord,
        lsst::afw::image::Exposure<float> const & exposure
    ) const {
        
        lsst::afw::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
        lsst::meas::base::CentroidResult result;
        result.x = center.getX() + _ctrl.param;
        result.y = center.getY() + _ctrl.param; 
        measRecord.set(_centroidKey, result); // better than NaN
    
        measRecord.set(_centroidKey, result);
        _flagHandler.setValue(measRecord, FAILURE, false);  // if we had a suspect flag, we'd set that instead
    }
    
    
    void fail(lsst::afw::table::SourceRecord & measRecord, lsst::meas::base::MeasurementError * error) const {
        _flagHandler.handleFailure(measRecord, error);
    }
    
    
    lsst::meas::base::CentroidResultKey _centroidKey;
    lsst::meas::base::FlagHandler _flagHandler;
    lsst::meas::base::SafeCentroidExtractor _centroidExtractor;
};

}}}
#endif // !LSST_MEAS_BASE_SillyCentroid_h_INCLUDED
