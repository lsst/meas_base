// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/table/arrays.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace base {

/**
 *  Configuration object for multiple-aperture flux algorithms
 */
class ApertureFluxControl {
public:

    ApertureFluxControl();

    LSST_CONTROL_FIELD(
        radii, std::vector<double>,
        "Radius (in pixels) of apertures."
    );

    LSST_CONTROL_FIELD(
        maxSincRadius, double,
        "Maximum radius (in pixels) for which the sinc algorithm should be used instead of the "
        "faster naive algorithm.  For elliptical apertures, this is the minor axis radius."
    );

    LSST_CONTROL_FIELD(
        shiftKernel, std::string,
        "Warping kernel used to shift Sinc photometry coefficients to different center positions"
    );

};

class ApertureFluxResult;

/**
 *  Base class for multiple-aperture photometry algorithms
 *
 *  ApertureFluxAlgorithm provides concrete implementations for single-aperture flux measurements,
 *  as well as a consistent interface and control object for its derived classes.  Currently, we only
 *  have one such derived class, CircularApertureFluxAlgorithm, but in the future we anticipate adding
 *  more derived classes for e.g. elliptical apertures, or apertures that are circular in sky coordinates
 *  rather than pixel coordinates.
 *
 *  Because the Result/ResultMapper system doesn't support the kind of flags ApertureFluxAlgorithm needs,
 *  we mostly don't use it here (though we use a Result object for the output of the single-aperture
 *  static methods), and that makes for a bit more boilerplate when creating Python plugins for
 *  these algorithms.  Hopefully that will be resolved in DM-1130.  In the meantime, this means that
 *  ApertureFluxAlgorithm and its subclasses cannot be wrapped using the WrappedSingleFramePlugin
 *  class, as most of the other C++ algorithm classes are.  Instead, manual plugin classes must be
 *  created that call the C++ measure() method.
 *
 *  Ultimately, ApertureFlux derived classes will fully replace the existing single-aperture plugins,
 *  SincFlux and NaiveFlux.  We can't do that yet, though, as we still can't use ApertureFlux outputs
 *  in slots (DM-1218).
 */
class ApertureFluxAlgorithm : public SimpleAlgorithm {
public:

    /// @copydoc PsfFluxAlgorithm::FlagBits

    enum {
        FAILURE=0,
        N_FLAGS
    };

    enum ResultFlagBits {
        APERTURE_TRUNCATED=0,
        SINC_COEFFS_TRUNCATED,
        RESULT_N_FLAGS
    };


    /// Typedef to the control object associated with this algorithm, defined above.
    typedef ApertureFluxControl Control;

    /// Result object returned by static methods.
    typedef ApertureFluxResult Result;

    //@{
    /**  Compute the flux (and optionally, uncertanties) within an aperture using Sinc photometry
     *
     *   The Sinc algorithm is slower than a naive aperture, but more accurate, in that it correctly
     *   handles sub-pixel aperture boundaries on well-sampled data.  This improved accuracy is most
     *   important for smaller apertures.
     *
     *   @param[in]   image                 Image or MaskedImage to be measured.  If a MaskedImage is
     *                                      provided, uncertainties will be returned as well as fluxes.
     *   @param[in]   ellipse               Ellipse that defines the outer boundary of the aperture.
     *   @param[in]   ctrl                  Control object.
     */
    template <typename T>
    static Result computeSincFlux(
        afw::image::Image<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    template <typename T>
    static Result computeSincFlux(
        afw::image::MaskedImage<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    //@}

    //@{
    /**  Compute the flux (and optionally, uncertanties) within an aperture using naive photometry
     *
     *   The naive algorithm just counts the flux in pixels whose centers lie within the aperture,
     *   ignoring the effects of sub-pixel aperture boundaries.
     *
     *   @param[in]   image                 Image or MaskedImage to be measured.  If a MaskedImage is
     *                                      provided, uncertainties will be returned as well as fluxes.
     *   @param[in]   ellipse               Ellipse that defines the outer boundary of the aperture.
     */
    template <typename T>
    static Result computeNaiveFlux(
        afw::image::Image<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    template <typename T>
    static Result computeNaiveFlux(
        afw::image::MaskedImage<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    //@}

    //@{
    /**  Compute the flux (and optionally, uncertanties) within an aperture using the algorithm
     *   determined by its size and the maxSincRadius control parameter.
     *
     *   This method delegates to computeSincFlux is the minor axis of the aperture is smaller than
     *   ctrl.maxSincRadius, and delegates to computeNaiveFlux otherwise.
     *
     *   @param[in]   image                 Image or MaskedImage to be measured.  If a MaskedImage is
     *                                      provided, uncertainties will be returned as well as fluxes.
     *   @param[in]   ellipse               Ellipse that defines the outer boundary of the aperture.
     *   @param[in]   ctrl                  Control object.
     */
    template <typename T>
    static Result computeFlux(
        afw::image::Image<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );

    template <typename T>
    static Result computeFlux(
        afw::image::MaskedImage<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    //@}

    /**
     *  Construct the algorithm and add its fields to the given Schema.
     */
    explicit ApertureFluxAlgorithm(
        Control const & ctrl,
        std::string const & name,
        afw::table::Schema & schema,
        daf::base::PropertySet & metadata
    );

    /**
     *  Measure the configured apertures on the given image.
     *
     *  Python plugins will delegate to this method.
     *
     *  @param[in,out] record      Record used to save outputs and retrieve positions.
     *  @param[in]     exposure    Image to be measured.
     */
    virtual void measure(
        afw::table::SourceRecord & record,
        afw::image::Exposure<float> const & exposure
    ) const = 0;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

    virtual ~ApertureFluxAlgorithm() {}

protected:

    struct FlagKeys {

        FlagKeys(std::string const & name, afw::table::Schema & schema, int index);

        afw::table::Key<afw::table::Flag> failed;
        afw::table::Key<afw::table::Flag> apertureTruncated;
        afw::table::Key<afw::table::Flag> sincCoeffsTruncated;
    };

    void copyResultToRecord(Result const & result, afw::table::SourceRecord & record, int index) const;

    Control const _ctrl;
    SafeCentroidExtractor _centroidExtractor;
    FlagHandler _flagHandler;
    afw::table::ArrayKey<Flux> _fluxKey;
    afw::table::ArrayKey<FluxErrElement> _fluxSigmaKey;
    std::vector<FlagKeys> _flagKeys;
};
struct ApertureFluxResult : public FluxResult {

    /// Return the flag value associated with the given bit
    bool getFlag(ApertureFluxAlgorithm::ResultFlagBits bit) const { return _flags[bit]; }

    /// Set the flag value associated with the given bit
    void setFlag(ApertureFluxAlgorithm::ResultFlagBits bit, bool value=true) { _flags[bit] = value; }

    /// Clear (i.e. set to false) the flag associated with the given bit
    void unsetFlag(ApertureFluxAlgorithm::ResultFlagBits bit) { _flags[bit] = false; }

private:

    std::bitset<ApertureFluxAlgorithm::RESULT_N_FLAGS> _flags;
};


}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_ApertureFlux_h_INCLUDED
