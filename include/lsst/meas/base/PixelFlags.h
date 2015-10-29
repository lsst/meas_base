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
#include <vector>

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle PixelFlagsAlgorithm's configuration
 */
class PixelFlagsControl {
public:
    LSST_CONTROL_FIELD(masksFpCenter, std::vector<std::string>,
                       "List of mask planes to be searched for which occur in the center of a footprint. "
                       "If any of the planes are found they will have a corresponding pixel flag set.");
    LSST_CONTROL_FIELD(masksFpAnywhere, std::vector<std::string>,
                       "List of mask planes to be searched for which occur anywhere within a footprint. "
                       "If any of the planes are found they will have a corresponding pixel flag set.");
    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    PixelFlagsControl() : masksFpCenter(), masksFpAnywhere() {}
};


/**
 *  @brief A measurement algorithm that gets mask bits from the exposure and sets flag bits to summarize which
 *         bits are set within a source's footprint.
 */
class PixelFlagsAlgorithm : public SimpleAlgorithm {
public:
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
        MeasurementError * error = NULL
    ) const;

    typedef std::map<std::string, afw::table::Key<afw::table::Flag>> KeyMap;

private:
    Control _ctrl;
    KeyMap _centerKeys;
    KeyMap _anyKeys;
    afw::table::Key<afw::table::Flag> _generalFailureKey;
    afw::table::Key<afw::table::Flag> _offImageKey;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_PixelFlags_h_INCLUDED
