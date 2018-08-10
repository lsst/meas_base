// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
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

#ifndef LSST_MEAS_BASE_ScaledApertureFlux_h_INCLUDED
#define LSST_MEAS_BASE_ScaledApertureFlux_h_INCLUDED

#include "lsst/afw/table.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/pex/config.h"

namespace lsst {
namespace meas {
namespace base {

class ScaledApertureFluxControl {
public:
    LSST_CONTROL_FIELD(
            shiftKernel, std::string,
            "Warping kernel used to shift Sinc photometry coefficients to different center positions");
    LSST_CONTROL_FIELD(scale, double, "Scaling factor of PSF FWHM for aperture radius.");

    // The default scaling factor is chosen such that scaled aperture
    // magnitudes are expected to be equal to Kron magnitudes, based on
    // measurements performed by Stephen Gwyn on WIRCam. See:
    // http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/wirwolf/docs/proc.html#photcal
    // http://www.cfht.hawaii.edu/fr/news/UM2013/presentations/Session10-SGwyn.pdf
    ScaledApertureFluxControl() : shiftKernel("lanczos5"), scale(3.14) {}
};

/**
 *  @brief Measure the instFlux in an aperture scaled to the PSF.
 *
 *  This algorithm performs a sinc aperture instFlux measurement where they size
 *  of the aperture is determined by multiplying the FWHM of the PSF by the
 *  scaling factor specified in the algorithm configuration.
 */
class ScaledApertureFluxAlgorithm : public SimpleAlgorithm {
public:
    typedef ScaledApertureFluxControl Control;
    typedef ApertureFluxResult Result;

    ScaledApertureFluxAlgorithm(Control const& control, std::string const& name, afw::table::Schema& schema);

    /**
     *  Measure the scaled aperture instFlux on the given image.
     *
     *  Python plugins will delegate to this method.
     *
     *  @param[in,out] record      Record used to save outputs and retrieve positions.
     *  @param[in]     exposure    Image to be measured.
     */
    virtual void measure(afw::table::SourceRecord& measRecord,
                         afw::image::Exposure<float> const& exposure) const override;

    virtual void fail(afw::table::SourceRecord& measRecord, MeasurementError* error = nullptr) const override;

private:
    Control _ctrl;
    FluxResultKey _instFluxResultKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
};

class ScaledApertureFluxTransform : public FluxTransform {
public:
    typedef ScaledApertureFluxControl Control;
    ScaledApertureFluxTransform(Control const& ctrl, std::string const& name,
                                afw::table::SchemaMapper& mapper);
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_ScaledApertureFlux_h_INCLUDED
