// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST.
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
#include "lsst/meas/base/Transform.h"

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


struct ApertureFluxResult;


/**
 *  Base class for multiple-aperture photometry algorithms
 *
 *  ApertureFluxAlgorithm serves as an intermediate base class for all aperture fluxes, which it assumes
 *  have that capability of measuring multiple apertures (even if they are not always configured to do
 *  so).
 *
 *  Concrete implementations for single-aperture flux measurements are provided as static methods,
 *  as well as a consistent interface and control object for its derived classes.  Currently, we only
 *  have one derived class, CircularApertureFluxAlgorithm, but in the future we anticipate adding
 *  more derived classes for e.g. elliptical apertures, or apertures that are circular in sky coordinates
 *  rather than pixel coordinates.
 */
class ApertureFluxAlgorithm : public SimpleAlgorithm {
public:

    /// @copydoc PsfFluxAlgorithm::FlagBits
    enum FlagBits {
        FAILURE=0,
        APERTURE_TRUNCATED,
        SINC_COEFFS_TRUNCATED,
        N_FLAGS
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

    virtual ~ApertureFluxAlgorithm() {}

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

    /// @copydoc BaseAlgorithm::fail
    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

protected:

    void copyResultToRecord(Result const & result, afw::table::SourceRecord & record, int index) const;

    FlagHandler const & getFlagHandler(int index) const { return _keys[index].flags; }

    Control const _ctrl;
    SafeCentroidExtractor _centroidExtractor;

private:

    struct Keys {
        FluxResultKey fluxKey;
        FlagHandler flags;

        Keys(afw::table::Schema & schema, std::string const & prefix, std::string const & doc, bool isSinc);
    };

    std::vector<Keys> _keys;
};


/**
 *  A Result struct for running an aperture flux algorithm with a single radius.
 *
 *  This simply extends FluxResult to add the appropriate error flags for aperture fluxes.
 */
struct ApertureFluxResult : public FluxResult {

    /// Return the flag value associated with the given bit
    bool getFlag(ApertureFluxAlgorithm::FlagBits bit) const { return _flags[bit]; }

    /// Set the flag value associated with the given bit
    void setFlag(ApertureFluxAlgorithm::FlagBits bit, bool value=true) { _flags[bit] = value; }

    /// Clear (i.e. set to false) the flag associated with the given bit
    void unsetFlag(ApertureFluxAlgorithm::FlagBits bit) { _flags[bit] = false; }

private:
    std::bitset<ApertureFluxAlgorithm::N_FLAGS> _flags;
};

/**
 *  Measurement transformation for aperture fluxes
 *
 *  Transforms fluxes with associated errors to magnitudes. Correctly handles
 *  multiple apertures. Flags are propagated from input to output.
 */
class ApertureFluxTransform : public BaseTransform {
public:
    typedef ApertureFluxControl Control;
    ApertureFluxTransform(Control const & ctrl, std::string const & name, afw::table::SchemaMapper & mapper);

    /*
     * @brief Perform transformation from inputCatalog to outputCatalog.
     *
     * @param[in]     inputCatalog   Source of data to be transformed
     * @param[in,out] outputCatalog  Container for transformed results
     * @param[in]     wcs            World coordinate system under which transformation will take place
     * @param[in]     calib          Photometric calibration under which transformation will take place
     * @throws        LengthError    Catalog sizes do not match
     */
    virtual void operator()(afw::table::SourceCatalog const & inputCatalog,
                            afw::table::BaseCatalog & outputCatalog,
                            afw::image::Wcs const & wcs,
                            afw::image::Calib const & calib) const;
private:
    std::vector<MagResultKey> _magKeys;
    Control _ctrl;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_ApertureFlux_h_INCLUDED
