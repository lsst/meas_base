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
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"
#include "lsst/meas/base/algorithms/SdssShapeImpl.h"

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
 *  @brief A measurement algorithm that gets mask bits from the exposure and sets flag bits to summarize.
 */
class PixelFlagsAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
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

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {{
                {"edge", "Could not use full PSF model image in fit because of proximity to exposure border"},
                {"interpolated", "Interpolated pixel in the source footprint"},
                {"interpolatedCenter", "Interpolated pixel in the source center"},
                {"saturated", "Saturated pixel in the source footprint"},
                {"saturatedCenter", "Saturated pixel in the source center"},
                {"cr", "Cosmic ray in the source footprint"},
                {"crCenter", "Cosmic ray in the source center"},
                {"bad", "Bad pixel in the source footprint"}

            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    typedef PixelFlagsControl Control;

    /**
     *  Result is the type returned by apply().
     */
    typedef Result0<PixelFlagsAlgorithm> Result;

    /**
     *  ResultMapper for PixelFlagsAlgorithm, which only defines a flag set out its output.
     */
    typedef ResultMapper0<PixelFlagsAlgorithm> ResultMapper;

    /**
     *  PixelFlagsAlgorithm only needs a footprint
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
     *  @brief Collect and summarize the masked image flags of the source footprint
     */
    template <typename T>
    static Result apply(
        afw::image::MaskedImage<T> const & image,
        afw::geom::Point2D const & position,
        afw::detection::Footprint const & footprint,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the PixelFlags to a single source using the Plugin API.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_PixelFlags_h_INCLUDED
