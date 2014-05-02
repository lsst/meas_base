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
 *  This algorithm is a flux measurement which simply sums the counts over the footprint
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle NaiveFluxAlgorithm's configuration
 */
class NaiveFluxControl {
public:
    LSST_CONTROL_FIELD(radius, double, "FIXME! NEVER DOCUMENTED!");


    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    NaiveFluxControl() : radius(7.0) {}
};


/**
 *  @brief A measurement algorithm that estimates flux using a linear least-squares fit with the Psf model
 *
 *  The NaiveFlux algorithm is extremely simple. It sums the pixels over the source footprint,
 *  within the radius specified in the config
 */
class NaiveFluxAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     *
     *  Note that we've included a final N_FLAGS value that isn't a valid flag; this is a common C++
     *  idiom for automatically counting the number of enum values, and it's required for Algorithms
     *  as the N_FLAGS value is used by the Result and ResultMapper objects.
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
    /// The control object contains the configuration parameters for this algorithm.
    typedef NaiveFluxControl Control;

    /**
     *  Result is the type returned by apply().  Because NaiveFluxAlgorithm only measures a flux and its
     *  uncertainty, we can use the single predefined component, FluxComponent, without any modification.
     *  Result1 is a template for algorithms with one result component, in addition to flags.
     */
    typedef Result1<NaiveFluxAlgorithm,FluxComponent> Result;

    /**
     *  The ResultMapper typedef here must exactly corresponds to the the Result typedef defined above:
     *  There is a ComponentMapper corresponding to each Component.
     */
    typedef ResultMapper1<NaiveFluxAlgorithm,FluxComponentMapper> ResultMapper;

    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  NaiveFluxAlgorithm only needs a centroid, so we use FootprintCentroidInput.
     */
    typedef FootprintCentroidInput Input; // type passed to apply in addition to Exposure.

    /**
     *  @brief Create an object that transfers Result values to a record associated with the given schema
     *
     *  This is called by the Plugin wrapper system to create a ResultMapper.  It's responsible for calling
     *  the ResultMapper constructor, forwarding the schema and prefix arguments and providing the correct
     *  values for the uncertainty arguments.
     */
    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Measure the flux of a source using the NaiveFlux algorithm.
     *
     *  This is the overload of apply() that does all the work, and it's designed to be as easy to use
     *  as possible outside the Plugin framework (since the Plugin framework calls the other one).  The
     *  arguments are all the things we need, and nothing more: we don't even pass a Footprint, since
     *  we wouldn't actually use it, and if we didn't need to get a Psf from the Exposure, we'd use
     *  MaskedImage instead.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the NaiveFlux to a single source using the Plugin API.
     *
     *  This is the version that will be called by both the SFM framework and the forced measurement
     *  framework, in single-object mode.  It will delegate to the other overload of apply().  Note that
     *  we can use the same implementation for both single-frame and forced measurement, because we require
     *  exactly the same inputs in both cases.  This is true for the vast majority of algorithms, but not
     *  all (extended source photometry is the notable exception).
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_NaiveFlux_h_INCLUDED
