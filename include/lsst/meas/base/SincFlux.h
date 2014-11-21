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

#ifndef LSST_MEAS_BASE_SincFlux_h_INCLUDED
#define LSST_MEAS_BASE_SincFlux_h_INCLUDED

/**
 *  @file lsst/meas/base/SincFlux.h
 *  This is the algorithm for SincFlux
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle SincFluxAlgorithm's configuration
 */
class SincFluxControl {
public:

    LSST_CONTROL_FIELD(radius1, double, "major axis of inner boundary (pixels)");
    LSST_CONTROL_FIELD(radius2, double, "major axis of outer boundary (pixels)");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    SincFluxControl() :
        radius1(0.0), radius2(7.0) {}
};


/**
 *  @brief An aperture photometry algorithm that uses sinc interpolation to handle fractional pixels.
 *
 *  While SincFluxAlgorithm supports elliptical apertures, the ellipticity and size of these apertures
 *  is fixed by config, so in practice this feature is not useful.
 *
 *  For larger apertures where the fractional nature of the aperture boundary is unimportant,
 *  NaiveFluxAlgorithm may be a better choice for performance reasons.
 *
 *  In the future, SincFluxAlgorithm and NaiveFluxAlgorithm will be merged into a single
 *  algorithm capable of computing multiple apertures (see DM-837) and handling elliptical
 *  apertures.
 */
class SincFluxAlgorithm : public SimpleAlgorithm {
public:

    enum FlagBits {
        FAILURE=FlagHandler::FAILURE,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef SincFluxControl Control;

    SincFluxAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

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
    FluxResultKey _fluxResultKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SincFlux_h_INCLUDED
