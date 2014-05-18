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

#ifndef LSST_MEAS_BASE_SdssShape_h_INCLUDED
#define LSST_MEAS_BASE_SdssShape_h_INCLUDED

/**
 *  @file lsst/meas/base/SdssShape.h
 *
 *  This file is one of two (the other is PsfFlux.h) intended to serve as an tutorial example on
 *  how to implement new Algorithms.  PsfFluxAlgorithm is a particularly simple algorithm, while
 *  SdssShapeAlgorithm is more complex.
 *
 *  See @ref measBaseImplementingNew for a general overview of the steps required.
 */

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle SdssShapeAlgorithm's configuration
 *
 *  @copydetails PsfFluxControl
 */
class SdssShapeControl {
public:
    LSST_CONTROL_FIELD(background, double, "Additional value to add to background");
    LSST_CONTROL_FIELD(maxIter, int, "Maximum number of iterations");
    LSST_CONTROL_FIELD(maxShift, double, "Maximum centroid shift");
    LSST_CONTROL_FIELD(tol1, float, "Convergence tolerance for e1,e2");
    LSST_CONTROL_FIELD(tol2, float, "Convergence tolerance for FWHM");

    /// @copydoc PsfFluxControl::PsfFluxControl
    SdssShapeControl() : background(0.0), maxIter(100), maxShift(), tol1(1E-5), tol2(1E-4) {}

};

/**
 *  @brief Additional results for SdssShapeAlgorithm
 *
 *  Unlike PsfFlux, some of SdssShape's outputs aren't handled by the standard FluxComponent,
 *  CentroidComponent, and ShapeComponent classes, so we have to define a Component class here
 *  to handle just those output which require special handling.
 *
 *  A corresponding ComponentMapper class is also required (see below).
 *
 *  Note: for what I guess are historical reasons, SdssShape computes covariance terms between the flux
 *  and the shape, but not between the flux and centroid or centroid and shape.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class SdssShapeExtras {
public:
    ShapeElement xy4;       ///< A fourth moment used in lensing (RHL needs to clarify; not in the old docs)
    ErrElement xy4Sigma;    ///< 1-Sigma uncertainty on xy4
    ErrElement flux_xx_Cov; ///< flux, xx term in the uncertainty covariance matrix
    ErrElement flux_yy_Cov; ///< flux, yy term in the uncertainty covariance matrix
    ErrElement flux_xy_Cov; ///< flux, xy term in the uncertainty covariance matrix

    SdssShapeExtras(); ///< Constructor; initializes everything to NaN
};

/**
 *  @brief Object that transfers additional SdssShapeAlgorithm results to afw::table records
 *
 *  Because we have custom outputs, we also have to define how to transfer those outputs to
 *  records.  We just follow the pattern established by the other ComponentMapper classes.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class SdssShapeExtrasMapper {
public:

    /**
     *  @brief Allocate fields in the schema and save keys for future use.
     *
     *  Unlike the standard ComponentMappers, SdssShapeExtrasMappers takes a Control instance
     *  as its third argument.  It doesn't actually need it, but this is a good pattern to
     *  establish as some algorithms' outputs *will* depend on the Control object's values, and
     *  the mapper needs to have some kind of third argument in order to work with
     *  @ref measBaseResultMapperTemplates.
     *
     *  All fields should start with the given prefix and an underscore.
     */
    SdssShapeExtrasMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        SdssShapeControl const & control=SdssShapeControl()
    );

    /// Transfer values from the result struct to the record.
    void apply(afw::table::BaseRecord & record, SdssShapeExtras const & result) const;

private:
    afw::table::Key<ShapeElement> _xy4;
    afw::table::Key<ErrElement> _xy4Sigma;
    afw::table::Key<ErrElement> _flux_xx_Cov;
    afw::table::Key<ErrElement> _flux_yy_Cov;
    afw::table::Key<ErrElement> _flux_xy_Cov;
};

/**
 *  @brief Measure the image moments of source using adaptive Gaussian weights.
 *
 *  This algorithm measures the weighted second moments of an image using a Gaussian weight function, which
 *  is iteratively updated to match the current weights.  If this iteration does not converge, it can fall
 *  back to using unweighted moments, which can be significantly noisier.
 *
 *  As an Algorithm class, all of SdssShapeAlgorithm's core functionality is available via static methods
 *  (in fact, there should be no reason to ever construct an instance).
 *
 *  Almost all of the implementation of SdssShapeAlgorithm is here and in SdssShapeAlgorithm.cc, but there
 *  are also a few key lines in the Swig .i file:
 *  @code
 *  %include "lsst/meas/base/SdssShape.h"
 *  %template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<float>;
 *  %template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<double>;
 *  %wrapMeasurementAlgorithm4(lsst::meas::base, SdssShapeAlgorithm, SdssShapeControl, FootprintCentroidInput,
 *                             ShapeComponent, CentroidComponent, FluxComponent, SdssShapeExtras)
 *  @endcode
 *  and in the pure Python layer:
 *  @code
 *  WrappedSingleFramePlugin.generate(SdssShapeAlgorithm)
 *  @endcode
 *  The former ensure the Algorithm class is fully wrapped via Swig (including @c %%template instantiations
 *  of its @c Result and @c ResultMapper classes), and the latter actually generates the Config class and
 *  the Plugin classes and registers them.
 */
class SdssShapeAlgorithm {
public:

    /// @copydoc PsfFluxAlgorithm::FlagBits
    enum FlagBits {
        UNWEIGHTED_BAD=0,
        UNWEIGHTED,
        SHIFT,
        MAXITER,
        N_FLAGS
    };

    /// @copydoc PsfFluxAlgorithm::getFlagDefinitions()
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {{
                {"unweightedBad", "Both weighted and unweighted moments were invalid"},
                {"unweighted", "Weighted moments converged to an invalid value; using unweighted moments"},
                {"shift", "centroid shifted by more than the maximum allowed amount"},
                {"maxIter", "Too many iterations in adaptive moments"}
            }};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef SdssShapeControl Control;

    /**
     *  This is the type returned by apply().  Because SdssShapeAlgorithm measures a flux, a centroid, a
     *  shape, and some other things, we concatenate all those together using the 4-parameter Result
     *  template. A Flags component is always included with these templates as well.
     */
    typedef Result4<
        SdssShapeAlgorithm,
        ShapeComponent,
        CentroidComponent,
        FluxComponent,
        SdssShapeExtras
    > Result;

    /// @copydoc PsfFluxAlgorithm::ResultMapper
    typedef ResultMapper4<
        SdssShapeAlgorithm,
        ShapeComponentMapper,
        CentroidComponentMapper,
        FluxComponentMapper,
        SdssShapeExtrasMapper
    > ResultMapper;

    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  SdssShapeAlgorithm needs a centroid, so we can use FootprintCentroidInput.
     */
    typedef FootprintCentroidInput Input;

    /// @copydoc PsfFluxAlgorithm::makeResultMapper
    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & name,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Measure the shape of a source using the SdssShape algorithm.
     *
     *  This is the overload of apply() that does all the work, and it's designed to be as easy to use
     *  as possible. The arguments contain everything needed, and nothing more: a MaskedImage is supplied
     *  (not an Exposure) because we just need the pixel data and mask.  Because the Plugin framework
     *  calls a different overload of apply() (the one that takes an Input object), it doesn't care what
     *  the signature of this one is.
     */
    template <typename T>
    static Result apply(
        afw::image::MaskedImage<T> const & image,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    /**
     *  @brief A limited implementation of SdssShape that runs on images with no variance plane, and
     *         hence reports no uncertainties.
     */
    template <typename T>
    static Result apply(
        afw::image::Image<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the SdssShape algorithm to a single source using the Plugin API.
     *
     *  This is the version that will be called by both the SFM framework.  It will delegate to the other
     *  overload of apply().
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssShape_h_INCLUDED
