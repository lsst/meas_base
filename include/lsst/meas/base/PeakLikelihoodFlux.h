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

#ifndef LSST_MEAS_BASE_PeakLikelihoodFlux_h_INCLUDED
#define LSST_MEAS_BASE_PeakLikelihoodFlux_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Transform.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief C++ control object for peak likelihood flux.
 *
 * Peak likelihood flux requires an image that has been filtered by convolving with its own PSF
 * (or an approximate model). It is equivalent to a PSF flux (e.g. PsfFluxAlgorithm) computed on
 * a non-preconvolved image.
 *
 * The PSF must be provided in the exposure, as it is used to compute a weighting factor.
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
 *  @brief A measurement algorithm that estimates the peak flux, using a filtered image
    which has been convolved with its own PSF.
 */
class PeakLikelihoodFluxAlgorithm : public SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    struct Flags;

    static std::size_t getFlagNumber(std::string const & name);
    static std::string const getFlagName(std::size_t flagNumber);

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef PeakLikelihoodFluxControl Control;

    PeakLikelihoodFluxAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

private:

    Control _ctrl;
    FluxResultKey _fluxResultKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
};

class PeakLikelihoodFluxTransform : public FluxTransform {
public:
    typedef PeakLikelihoodFluxControl Control;
    PeakLikelihoodFluxTransform(Control const & ctrl, std::string const & name,
                                afw::table::SchemaMapper & mapper) : FluxTransform{name, mapper} { }
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_PeakLikelihoodFlux_h_INCLUDED
