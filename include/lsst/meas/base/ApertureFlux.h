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
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

class ApertureFluxControl {
public:

    ApertureFluxControl() : maxSincRadius(10.0), shiftKernel("lanczos5") {}

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

/**
 *  @brief A reusable component for result structs that contain an array of flux measurements.
 *
 *  Flux measurements and their errors should always be in DN.
 */
struct ApertureFluxComponent {
    ndarray::Array<Flux,1,1> flux; ///< Measured flux in DN.
    ndarray::Array<FluxErrElement,1,1> fluxSigma; ///< 1-Sigma error (sqrt of variance) on flux in DN.

    /// Return a FluxComponent corresponding to the ith flux
    FluxComponent const get(int i) const { return FluxComponent(flux[i], fluxSigma[i]); }

    /// Return set the ith flux and its uncertainty to the values in the given FluxComponent
    void set(int i, FluxComponent const & c) const {
        flux[i] = c.flux;
        fluxSigma[i] = c.fluxSigma;
    }

    /// Default constructor; arrays will remain empty
    ApertureFluxComponent() {}

    /// Constructor; arrays will be allocated to the given size and filled with NaN
    explicit ApertureFluxComponent(int size);
};

/**
 *  @brief An object that transfers values from ApertureFluxComponent to afw::table::BaseRecord
 *
 *  This should be included in one of @ref measBaseResultMapperTemplates to correspond with using
 *  ApertureFluxComponent in the same position in one of @ref measBaseResultTemplates, and will otherwise
 *  not be used directly by users.
 */
class ApertureFluxComponentMapper {
public:

    /**
     *  @brief Construct the mapper, adding fields to the given schema and saving their keys
     *
     *  The given prefix will form the first part of all fields, and the uncertainty argument
     *  sets which uncertainty fields will be added to the schema and transferred during apply().
     */
    ApertureFluxComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        std::vector<double> const & radii
    );

    /// Transfer values from the result struct to the record
    void apply(afw::table::BaseRecord & record, ApertureFluxComponent const & result) const;

private:
    afw::table::ArrayKey<Flux> _flux;
    afw::table::ArrayKey<FluxErrElement> _fluxSigma;
};

class ApertureFluxAlgorithm {
public:

    typedef ApertureFluxControl Control;

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
    static Flux computeSincFlux(
        afw::image::Image<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    template <typename T>
    static FluxComponent computeSincFlux(
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
    static Flux computeNaiveFlux(
        afw::image::Image<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    );
    template <typename T>
    static FluxComponent computeNaiveFlux(
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
    static Flux computeFlux(
        afw::image::Image<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    ) {
        return (afw::geom::ellipses::Axes(ellipse.getCore()).getB() <= ctrl.maxSincRadius)
            ? computeSincRadius(image, ellipse, ctrl)
            : computeNaiveRadius(image, ellipse, ctrl);
    }
    template <typename T>
    static FluxComponent computeFlux(
        afw::image::MaskedImage<T> const & image,
        afw::geom::ellipses::Ellipse const & ellipse,
        Control const & ctrl=Control()
    ) {
        return (afw::geom::ellipses::Axes(ellipse.getCore()).getB() <= ctrl.maxSincRadius)
            ? computeSincRadius(image, ellipse, ctrl)
            : computeNaiveRadius(image, ellipse, ctrl);
    }
    //@}

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_ApertureFlux_h_INCLUDED
