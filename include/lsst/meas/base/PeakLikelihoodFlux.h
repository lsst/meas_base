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

#ifndef LSST_MEAS_BASE_PeakLikelihoodFlux_h_INCLUDED
#define LSST_MEAS_BASE_PeakLikelihoodFlux_h_INCLUDED

/**
 *  @file lsst/meas/base/PeakLikelihoodFlux.h
 *
 *  This algorithm is a flux measurement which simply sums the counts over the footprint
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief C++ control object for peak likelihood flux.
 *
 * Peak likelihood flux requires an image that has been filtered by convolving with its own PSF
 * (or an approximate model). The PSF must be provided in the exposure, as it is used to compute
 * a weighting factor.
 * 
 * Flux and error are computed as follows:
 * * flux = sum(unfiltered image * PSF) / sum(PSF^2)
 *        = value of peak of filtered source / sum(PSF^2)
 * * err  = sqrt(sum(unfiltered variance * PSF^2) / sum(PSF^2)^2)
 *        = sqrt(value of filtered variance at peak / sum(PSF^2)^2)
 * * The pixels in the image are samples of a band-limited function, and by using
 *   a sinc interpolation (via a warping kernel) we can evaluate this function at any point.
 *   We use this technique to compute the peak of the function, which is assumed to be
 *   at the centroid of the filtered source.
 */
class PeakLikelihoodFluxControl {
public:

    LSST_CONTROL_FIELD(warpingKernelName, std::string,
        "Name of warping kernel (e.g. \"lanczos4\") used to compute the peak");

    PeakLikelihoodFluxControl() : warpingKernelName("lanczos4") {}
};
/**
 *  @brief A measurement algorithm that estimates flux using a linear least-squares fit with the Psf model
 *
 *  The PeakLikelihoodFlux algorithm is extremely simple. It sums the pixels over the source footprint,
 *  within the radius specified in the config
 */
class PeakLikelihoodFluxAlgorithm {
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
    typedef PeakLikelihoodFluxControl Control;

    /**
     *  Result is the type returned by apply().  Because PeakLikelihoodFluxAlgorithm only measures a flux and its
     *  uncertainty, we can use the single predefined component, FluxComponent, without any modification.
     */
    typedef Result1<PeakLikelihoodFluxAlgorithm,FluxComponent> Result;

    /**
     *  The ResultMapper typedef here must exactly corresponds to the the Result typedef defined above:
     *  There is a ComponentMapper corresponding to each Component.
     */
    typedef ResultMapper1<PeakLikelihoodFluxAlgorithm,FluxComponentMapper> ResultMapper;

    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  PeakLikelihoodFluxAlgorithm only needs a centroid, so we use FootprintCentroidInput.
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
     *  @brief Measure the flux of a source using the PeakLikelihoodFlux algorithm.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the PeakLikelihoodFlux to a single source using the Plugin API.
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

#endif // !LSST_MEAS_BASE_PeakLikelihoodFlux_h_INCLUDED
