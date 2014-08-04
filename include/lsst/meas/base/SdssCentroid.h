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

#ifndef LSST_MEAS_BASE_SdssCentroid_h_INCLUDED
#define LSST_MEAS_BASE_SdssCentroid_h_INCLUDED

/**
 *  @file lsst/meas/base/SdssCentroid.h
 *
 *  This is implementation of SdssCentroid using the meas_base measurement framework
 *
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle SdssCentroidAlgorithm's configuration
 *
 */
class SdssCentroidControl {
public:

    LSST_CONTROL_FIELD(binmax, int, "maximum allowed binning");
    LSST_CONTROL_FIELD(peakMin, double, "if the peak's less than this insist on binning at least once");
    LSST_CONTROL_FIELD(wfac, double, "fiddle factor for adjusting the binning");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */

    SdssCentroidControl() : binmax(16), peakMin(-1.0), wfac(1.5) {}
};

/**
 *  @brief The Sdss Centroid Algorithm
 */
class SdssCentroidAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
        NO_PSF,
        EDGE,
        BAD_DATA,
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {{
                {"noPsf", "Exposure has no attached Psf"},
                {"badData", "Algorithm could not measure this data"},
                {"edge", "Object too close to edge"}
            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef SdssCentroidControl Control;

    /**
     *  This is the type returned by apply().  SdssCentroidAlgorithm returns only a Point and Errors
     */
    typedef Result1<
        SdssCentroidAlgorithm,
        CentroidComponent 
    > Result;

    /// @copydoc PsfFluxAlgorithm::ResultMapper
    typedef ResultMapper1<
        SdssCentroidAlgorithm,
        CentroidComponentMapper
    > ResultMapper;


    /**
     *  Input from the measurement framework to the algorithm.
     *  SdssCentroidAlgorithm needs only a beginning centroid and footprint,
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
     *  @brief Measure the centroid of a source using the SdssCentroid algorithm.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the SdssCentroid to a single source using the Plugin API.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssCentroid_h_INCLUDED
