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

#ifndef LSST_MEAS_BASE_GaussianFlux_h_INCLUDED
#define LSST_MEAS_BASE_GaussianFlux_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle GaussianFluxAlgorithm's configuration
 */
class GaussianFluxControl {
public:
    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    GaussianFluxControl() :
        background(0.0)
    {}
};

/**
 *  @brief A measurement algorithm that estimates flux using an elliptical Gaussian weight.
 *
 *  This algorithm computes flux as the dot product of an elliptical Gaussian weight function
 *  with the image.  The size and ellipticity of the weight function is determined using the
 *  SdssShape algorithm, or retreived from a named field.
 */
class GaussianFluxAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
        N_FLAGS=0
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
    typedef GaussianFluxControl Control;

    /**
     *  Result is the type returned by apply().  Because GaussianFluxAlgorithm only measures a flux and its
     *  uncertainty, we can use the single predefined component, FluxComponent, without modification.
     */
    typedef Result1<GaussianFluxAlgorithm,FluxComponent> Result;

    /**
     *  Use the FluxComponentMapper to map algorithm Result to output catalog
     */
    typedef ResultMapper1<GaussianFluxAlgorithm,FluxComponentMapper> ResultMapper;

    /**
     *  GaussianFluxAlgorithm only needs a centroid and footprint as input.
     */
    typedef FootprintCentroidShapeInput Input; // type passed to apply in addition to Exposure.

    /**
     *  @brief Create an object that transfers Result values to a record associated with the given schema
     */
    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Measure the flux of a source using the GaussianFlux algorithm.
     */
    template <typename T>
    static void apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & centroid,
        afw::geom::ellipses::Quadrupole const & shape,
        Result & result,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the GaussianFlux to a single source using the Plugin API.
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

#endif // !LSST_MEAS_BASE_GaussianFlux_h_INCLUDED
