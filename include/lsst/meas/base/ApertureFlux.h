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

#ifndef LSST_MEAS_BASE_ApertureFlux_h_INCLUDED
#define LSST_MEAS_BASE_ApertureFlux_h_INCLUDED

/**
 *  @file lsst/meas/base/ApertureFlux.h
 *
 */
 
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"
// #include "lsst/meas/base/ApFluxComponent.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle ApertureFluxAlgorithm's configuration
 */
class ApertureFluxControl {
public:
    LSST_CONTROL_FIELD(radii, std::vector<double>, "vector of radii for apertures (in pixels)");

    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be excluded from the fit");
    /**
    *  @brief Default constructor
    *
    *  All control classes should define a default constructor that sets all fields to their default values.
    */
    ApertureFluxControl();
};
 
/**
 *  @brief An object that transfers values from FluxComponent to afw::table::BaseRecord
 *
 *  This should be included in one of @ref measBaseResultMapperTemplates to correspond with using
 *  FluxComponent in the same position in one of @ref measBaseResultTemplates, and will otherwise
 *  not be used directly by users.
 */
class ApFluxComponentMapper {
public:
 
    /**
     *  @brief Construct the mapper, adding fields to the given schema and saving their keys
     *  sets which uncertainty fields will be added to the schema and transferred during apply().
     */
    ApFluxComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        UncertaintyEnum uncertainty
    );

    /// Transfer values from the result struct to the record
    void apply(afw::table::BaseRecord & record, FluxComponent const & result) const;

private:
    afw::table::Key<Flux> _flux;
    afw::table::Key<ErrElement> _fluxSigma;
};
 
 /**
  *  @brief Additional results for ApertureFluxAlgorithm
  *
  *  ApertureFlux keeps an a vector of FluxComponents of the same length as the config.radii
  *  parameter.  So it can use the standard FluxComponent, but not the standard mapper.
  */
 
class ApertureFluxExtras {
public:
    std::vector< boost::shared_ptr<lsst::meas::base::FluxComponent> > fluxComponentVector; //changed pgee
    ApertureFluxExtras(); ///< Constructor; initializes everything to NaN
};


/**
 *  @brief Object that transfers additional ApertureFlux flux component vector to an output
 *  source catalog. 
 */
class ApertureFluxExtrasMapper {
public:

    /**
     *  @brief Allocate fields in the schema and save keys for future use.
     *  The output fluxes are named with a number postpended to the field name, equal to
     *  the index of the radius in config.radii and to the component in fluxComponentVector.
     */
    ApertureFluxExtrasMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ApertureFluxControl const & control=ApertureFluxControl() 
    );

    /// Transfer values from the result struct to the record.
    void apply(afw::table::BaseRecord & record, ApertureFluxExtras const & result) const;

private:
    afw::table::Key<int> _nApertures;
    std::vector<FluxComponentMapper> _fluxComponentMapperVector; //changed pgee
};

/**
 *  @brief A measurement algorithm that measures the flux within each aperture radius
 *  in the config.radii list, and outputs that number of FluxComponent results to the catalog.
 */
class ApertureFluxAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
        NO_PSF=0,
        NO_GOOD_PIXELS,
        EDGE,
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {{
                {"noPsf", "No Psf object attached to the Exposure object being measured"},
                {"noGoodPixels", "No usable pixels in fit region"},
                {"edge", "At least one of the aperture extends to or past the edge"}
            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef ApertureFluxControl Control;

    /**
     *  This is the type returned by apply().  Because ApertureFluxAlgorithm measure multiple fluxes,
     *  the Result object must be a custom vector of FluxComponents
     */
    typedef Result1<
        ApertureFluxAlgorithm,
        ApertureFluxExtras
    > Result;

    /// @copydoc PsfFluxAlgorithm::ResultMapper
    typedef ResultMapper1<
        ApertureFluxAlgorithm,
        ApertureFluxExtrasMapper
    > ResultMapper;


    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  ApertureFluxAlgorithm only needs a centroid, so we use FootprintCentroidInput.
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
     *  @brief Measure multiple fluxes of a source using the ApertureFlux algorithm.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the ApertureFlux to a single source using the Plugin API.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_ApertureFlux_h_INCLUDED
