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

#ifndef LSST_MEAS_PLUGINS_SillyCentroid_h_INCLUDED
#define LSST_MEAS_PLUGINS_SillyCentroid_h_INCLUDED

/**
 *  @file lsst/meas/base/SillyCentroid.h
 *
 *  This implements the SillyCentroid algorithm within the meas_base measurement framework
 *
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace test { namespace foo { namespace bar {

/**
 *  @brief A C++ control class to handle SillyCentroidAlgorithm's configuration
 *
 */
class SillyCentroidControl {
public:

    LSST_CONTROL_FIELD(param, int, "Difference to apply to y");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */

    SillyCentroidControl() : param(0) {}
};

/**
 *  @brief The Silly Centroid Algorithm
 */
class SillyCentroidAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
        EDGE,
        BAD_DATA,
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<lsst::meas::base::FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<lsst::meas::base::FlagDef,N_FLAGS> const flagDefs = {{
                {"badData", "Algorithm could not measure this data"},
                {"edge", "Object too close to edge"}
            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef SillyCentroidControl Control;

    /**
     *  This is the type returned by apply().  SillyCentroidAlgorithm returns only a Point and Errors
     */
    typedef lsst::meas::base::Result1<
        SillyCentroidAlgorithm,
        lsst::meas::base::CentroidComponent 
    > Result;

    /// @copydoc PsfFluxAlgorithm::ResultMapper
    typedef lsst::meas::base::ResultMapper1<
        SillyCentroidAlgorithm,
        lsst::meas::base::CentroidComponentMapper
    > ResultMapper;


    /**
     *  Input from the measurement framework to the algorithm.
     *  SillyCentroidAlgorithm needs only a beginning centroid and footprint.
     */
    typedef lsst::meas::base::FootprintCentroidInput Input; // type passed to apply in addition to Exposure.

    /**
     *  @brief Create an object that transfers Result values to a record associated with the given schema
     */
    static ResultMapper makeResultMapper(
        lsst::afw::table::Schema & schema,
        std::string const & prefix,
        Control const & ctrl=Control()
    ) {
        return ResultMapper(schema, prefix, lsst::meas::base::FULL_COVARIANCE);
    }

    /**
     *  @brief Measure the centroid of a source using the SillyCentroid algorithm.
     */
    static void apply(
        lsst::afw::image::Exposure<double> const & exposure,
        lsst::afw::geom::Point2D const & position,
        Result & result,
        Control const & ctrl=Control()
    ) {
        lsst::afw::geom::Point2D newpos = position + lsst::afw::geom::Extent2D(ctrl.param, ctrl.param);
        result.x = newpos.getX();
        result.y = newpos.getY();
    }
    /**
     *  @brief Apply the SillyCentroid to a single source using the Plugin API.
     */
    static void apply(
        lsst::afw::image::Exposure<double> const & exposure,
        Input const & inputs,
        Result & result,
        Control const & ctrl=Control()
    ) { 
        apply(exposure, inputs.position, result, ctrl);
    }

    static void apply(
        lsst::afw::image::Exposure<float> const & exposure,
        lsst::afw::geom::Point2D const & position,
        Result & result,
        Control const & ctrl=Control()
    ) {
        lsst::afw::geom::Point2D newpos = position + lsst::afw::geom::Extent2D(ctrl.param, ctrl.param);
        result.x = newpos.getX();
        result.y = newpos.getY();
    }
    /**
     *  @brief Apply the SillyCentroid to a single source using the Plugin API.
     */
    static void apply(
        lsst::afw::image::Exposure<float> const & exposure,
        Input const & inputs,
        Result & result,
        Control const & ctrl=Control()
    ) { 
        apply(exposure, inputs.position, result, ctrl);
    }
};

}}} // namespace lsst::meas::base

#endif
