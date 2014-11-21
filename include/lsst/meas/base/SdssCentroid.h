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

#ifndef LSST_MEAS_BASE_SdssCentroid_h_INCLUDED
#define LSST_MEAS_BASE_SdssCentroid_h_INCLUDED

/**
 *  @file lsst/meas/base/SdssCentroid.h
 *
 *  This implements the SdssCentroid algorithm within the meas_base measurement framework
 *
 */

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle SdssCentroidAlgorithm's configuration
 *
 */
class SdssCentroidControl {
public:

    LSST_CONTROL_FIELD(binmax, int, "maximum allowed binning");
    LSST_CONTROL_FIELD(peakMin, double, "if the peak's less than this insist on binning at least once");
    LSST_CONTROL_FIELD(wfac, double, "fiddle factor for adjusting the binning");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */

    SdssCentroidControl() : binmax(16), peakMin(-1.0), wfac(1.5) {}
};

/**
 *  @brief The Sdss Centroid Algorithm
 */
class SdssCentroidAlgorithm : public SimpleAlgorithm {
public:

    enum FlagBits {
        FAILURE=FlagHandler::FAILURE,
        EDGE,
        BAD_DATA,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef SdssCentroidControl Control;

    SdssCentroidAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

private:

    // These are private so they doesn't shadow the other overloads in base classes;
    // we can still call it via the public method on the base class.  We could have
    // used a using declaration instead, but Swig had trouble with that here.

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

    Control _ctrl;
    CentroidResultKey _centroidKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssCentroid_h_INCLUDED
