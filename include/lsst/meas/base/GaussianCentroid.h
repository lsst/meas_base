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

#ifndef LSST_MEAS_BASE_GaussianCentroid_h_INCLUDED
#define LSST_MEAS_BASE_GaussianCentroid_h_INCLUDED

#include <cmath>
#include <vector>
#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

struct FittedModel {
    enum { PEAK = 0, SKY, X0, Y0, SIGMA, NPARAM };
    
    enum {
        BAD_GUESS = -11,
        TOO_FEW = -12,
        CHI_SQUARED = -13,
        RANGE = -14,
        BAD_WIDTH = -15,
        LOST = -16,
        DIAGONAL = -17,
        BAD_A = -18,
        CONVERGE = 1,
        ITERATE = 2,
        ALMOST = 3,
        POOR = 4
    };

    FittedModel(int status_, std::vector<double> params_, int iter_=0, double flamd_=0, double chnew_=0) :
        status(status_), params(params_), iter(iter_), flamd(flamd_), chnew(chnew_) { }
    int status;
    std::vector<double> params;
    int iter;
    double flamd;
    double chnew;
};


/**
 *  @brief A C++ control class to handle GaussianCentroidAlgorithm's configuration
 *
 *  At present, GaussianCentroidAlgorithm has no configuration options.
 */
class GaussianCentroidControl {
public:
    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    GaussianCentroidControl() {}
};

/**
 *  @brief A class that calculates a centroid by fitting a circular Gaussian to the image.
 */
class GaussianCentroidAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum FlagBits {
        EDGE,
        NO_PEAK,
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {{
                {"edge", "Near edge of image"},
                {"noPeak", "Fitted Centroid has a negative peak"}
            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef GaussianCentroidControl Control;

    /**
     *  This is the type returned by apply().
     */
    typedef Result1<
        GaussianCentroidAlgorithm,
        CentroidComponent
    > Result;

    /// @copydoc PsfFluxAlgorithm::ResultMapper
    typedef ResultMapper1<
        GaussianCentroidAlgorithm,
        CentroidComponentMapper
    > ResultMapper;


    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  GaussianCentroidAlgorithm only needs a centroid, so we use
     *  FootprintCentroidInput.
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
     *  @brief Measure the centroid of a source using the GaussianCentroid algorithm.
     */
    template <typename T>
    static void apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & position,
        Result & result,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the GaussianCentroid to a single source using the Plugin API.
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

#endif // !LSST_MEAS_BASE_GaussianCentroid_h_INCLUDED
