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

#ifndef LSST_MEAS_BASE_NaiveFlux_h_INCLUDED
#define LSST_MEAS_BASE_NaiveFlux_h_INCLUDED

/**
 *  @file lsst/meas/base/NaiveFlux.h
 *
 *  This algorithm estimates flux by summing the pixel values at the given center
 *  over a circular aperture.  The pixel radius is provided by the configuration.
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle NaiveFluxAlgorithm's configuration
 */
class NaiveFluxControl {
public:
    LSST_CONTROL_FIELD(radius, double, "Size of the circular aperture, in pixels");


    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    NaiveFluxControl() : radius(7.0) {}
};


/**
 *  @brief A measurement algorithm that sums pixel values over a circular aperture.
 *
 *  The NaiveFlux algorithm is extremely simple. It sums the pixels over the source footprint,
 *  within the radius specified in the config.  This does not account for the fact that the
 *  aperture boundary involves fractional pixels, as SincFluxAlgorithm, but it is much faster
 *  and hence should be used for large apertures where the fractional boundary pixels are
 *  less important.
 *
 *  In the future, SincFluxAlgorithm and NaiveFluxAlgorithm will be merged into a single
 *  algorithm capable of computing multiple apertures (see DM-837) and handling elliptical
 *  apertures.
 */
class NaiveFluxAlgorithm : public SimpleAlgorithm {
public:

    enum FlagBits {
        FAILURE=FlagHandler::FAILURE,
        EDGE,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef NaiveFluxControl Control;

    NaiveFluxAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

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

#endif // !LSST_MEAS_BASE_NaiveFlux_h_INCLUDED
