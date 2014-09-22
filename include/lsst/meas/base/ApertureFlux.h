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
