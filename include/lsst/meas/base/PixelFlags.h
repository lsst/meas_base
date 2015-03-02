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

#ifndef LSST_MEAS_BASE_PixelFlags_h_INCLUDED
#define LSST_MEAS_BASE_PixelFlags_h_INCLUDED

/**
 *  @file lsst/meas/base/PixelFlags.h
 *  This is the algorithm for PixelFlags
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle PixelFlagsAlgorithm's configuration
 */
class PixelFlagsControl {
public:

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    PixelFlagsControl() {}
};


/**
 *  @brief A measurement algorithm that gets mask bits from the exposure and sets flag bits to summarize which
 *         bits are set within a source's footprint.
 */
class PixelFlagsAlgorithm : public SimpleAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum {
        FAILURE=FlagHandler::FAILURE,
        EDGE,
        INTERPOLATED,
        INTERPOLATED_CENTER,
        SATURATED,
        SATURATED_CENTER,
        CR,
        CR_CENTER,
        BAD,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef PixelFlagsControl Control;

    PixelFlagsAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

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

#endif // !LSST_MEAS_BASE_PixelFlags_h_INCLUDED
