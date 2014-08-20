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
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"
#include "lsst/meas/base/algorithms/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle SincFluxAlgorithm's configuration
 */
class SincFluxControl {
public:

    LSST_CONTROL_FIELD(radius1, double, "major axis of inner boundary (pixels)");
    LSST_CONTROL_FIELD(radius2, double, "major axis of outer boundary (pixels)");
    LSST_CONTROL_FIELD(angle, double, "measured from x anti-clockwise; radians");
    LSST_CONTROL_FIELD(ellipticity, double, "1 - b/a");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    SincFluxControl() :
        radius1(0.0), radius2(7.0), angle(0.0), ellipticity(0.0) {}
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
class SincFluxAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    typedef SincFluxControl Control;

    /**
     *  Result is the type returned by apply().  Because SincFluxAlgorithm only measures a flux and its
     *  uncertainty, we can use the single predefined component, FluxComponent, without any modification.
     */
    typedef Result1<SincFluxAlgorithm,FluxComponent> Result;

    /**
     *  The ResultMapper typedef here must exactly corresponds to the the Result typedef defined above:
     *  There is a ComponentMapper corresponding to each Component.
     */
    typedef ResultMapper1<SincFluxAlgorithm,FluxComponentMapper> ResultMapper;

    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  SincFluxAlgorithm only needs a centroid, so we use FootprintCentroidInput.
     */
    typedef FootprintCentroidInput Input; // type passed to apply in addition to Exposure.

    /**
     *  @brief Create an object that transfers Result values to a record associated with the given schema
     */
    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Measure the flux of a source using the SincFlux algorithm.
     */
    template <typename T>
    static void apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & position,
        Result & result,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the SincFlux to a single source using the Plugin API.
     */
    template <typename T>
    static void apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Result & result,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SincFlux_h_INCLUDED
