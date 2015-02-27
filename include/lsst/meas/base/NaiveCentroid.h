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

#ifndef LSST_MEAS_BASE_NaiveCentroid_h_INCLUDED
#define LSST_MEAS_BASE_NaiveCentroid_h_INCLUDED

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace base {
/**
 *  @brief A C++ control class to handle NaiveCentroidAlgorithm's configuration
 */
class NaiveCentroidControl {
public:
    LSST_CONTROL_FIELD(background, double, "Value to subtract from the image pixel values");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    NaiveCentroidControl() : background(0.0) {}
};

/**
 *  @brief A class that calculates a centroid as a simple unweighted first moment
 *         of the 3x3 region around a pixel.
 *
 *   A fixed background (set via config) may optionally be subtracted.
 *   This algorithm does not currently report an error, but it probably should.
 */

class NaiveCentroidAlgorithm : public SimpleAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum {
        FAILURE=FlagHandler::FAILURE,
        NO_COUNTS,
        EDGE,
        N_FLAGS
    };

    typedef NaiveCentroidControl Control;

    NaiveCentroidAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

private:

    Control _ctrl;
    CentroidResultKey _centroidKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_NaiveCentroid_h_INCLUDED
