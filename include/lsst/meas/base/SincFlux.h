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
 *
 *  This file is one of two (the other is SdssShape.h) intended to serve as an tutorial example on
 *  how to implement new Algorithms.  SincFluxAlgorithm is a particularly simple algorithm, while
 *  SdssShapeAlgorithm is more complex.
 *
 *  See @ref measBaseImplementingNew for a general overview of the steps required.
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"
#include "lsst/meas/base/algorithms/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle SincFluxAlgorithm's configuration
 *
 *  In C++, we define Control objects to handle configuration information.  Using the LSST_CONTROL_FIELD
 *  macro and lsst.pex.config.wrap.makeConfigClass, we can turn these into more full-featured Config classes
 *  in Python.  While the user will usually interact with the Config class, the plugin wrapper system will
 *  turn Config instances into Control instances when passing them to C++.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
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
 *  @brief A measurement algorithm that estimates flux using a linear least-squares fit with the Psf model
 *
 *  The SincFlux algorithm is extremely simple: we do a least-squares fit of the Psf model (evaluated
 *  at a given position) to the data.  For point sources, this provides the optimal flux measurement
 *  in the limit where the Psf model is correct.  We do not use per-pixel weights in the fit by default
 *  (see SincFluxControl::usePixelWeights), as this results in bright stars being fit with a different
 *  effective profile than faint stairs.
 *
 *  As one of the simplest Algorithms, SincFlux is documented to serve as an example in implementing new
 *  algorithms.  For an overview of the interface Algorithms should adhere to, see
 *  @ref measBaseAlgorithmConcept.
 *
 *  As an Algorithm class, all of SincFluxAlgorithm's core functionality is available via static methods
 *  (in fact, there should be no reason to ever construct an instance).
 *
 *  Almost all of the implementation of SincFluxAlgorithm is here and in SincFluxAlgorithm.cc, but there
 *  are also a few key lines in the Swig .i file:
 *  @code
 *  %include "lsst/meas/base/SincFlux.h"
 *  %template(apply) lsst::meas::base::SincFluxAlgorithm::apply<float>;
 *  %template(apply) lsst::meas::base::SincFluxAlgorithm::apply<double>;
 *  %wrapMeasurementAlgorithm1(lsst::meas::base, SincFluxAlgorithm, SincFluxControl, FootprintCentroidInput,
 *                             FluxComponent)
 *  @endcode
 *  and in the pure Python layer:
 *  @code
 *  WrappedSingleFramePlugin.generate(SincFluxAlgorithm)
 *  @endcode
 *  The former ensure the Algorithm class is fully wrapped via Swig (including @c %%template instantiations
 *  of its @c Result and @c ResultMapper classes), and the latter actually generates the Config class and
 *  the Plugin classes and registers them.
 */
class SincFluxAlgorithm {
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
        NO_PSF=0,
        NO_GOOD_PIXELS,
        EDGE,
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     *
     *  Each element of the returned array should correspond to one of the FlagBits enum values, but the
     *  names should follow conventions; FlagBits should be ALL_CAPS_WITH_UNDERSCORES, while FlagDef names
     *  should be camelCaseStartingWithLowercase.  @sa FlagsComponentMapper.
     *
     *  The implementation of getFlagDefinitions() should generally go in the header file so it is easy
     *  to keep in sync with the FlagBits enum.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {{
                {"noPsf", "No Psf object attached to the Exposure object being measured"},
                {"noGoodPixels", "No usable pixels in fit region"},
                {"edge", "Could not use full PSF model image in fit because of proximity to exposure border"}
            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef SincFluxControl Control;

    /**
     *  Result is the type returned by apply().  Because SincFluxAlgorithm only measures a flux and its
     *  uncertainty, we can use the single predefined component, FluxComponent, without any modification.
     *  Result1 is a template for algorithms with one result component, in addition to flags.
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
     *  @brief Measure the flux of a source using the SincFlux algorithm.
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
     *  @brief Apply the SincFlux to a single source using the Plugin API.
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

#endif // !LSST_MEAS_BASE_SincFlux_h_INCLUDED
