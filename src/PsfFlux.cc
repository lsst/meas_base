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

#include "boost/array.hpp"

#include "ndarray/eigen.h"

#include "lsst/afw/table/Source.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/meas/base/PsfFlux.h"

namespace lsst { namespace meas { namespace base {

PsfFluxAlgorithm::PsfFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxResultKey(
        FluxResultKey::addFields(schema, name, "flux derived from linear least-squares fit of PSF model")
    )
{
    static boost::array<FlagDefinition,N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag"},
        {"flag_noGoodPixels", "not enough non-rejected pixels in data to attempt the fit"},
        {"flag_edge", "object was too close to the edge of the image to use the full PSF model"}
    }};
    _flagHandler = FlagHandler(schema, name, flagDefs.begin(), flagDefs.end());
}

void PsfFluxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            FatalAlgorithmError,
            "PsfFlux algorithm requires a Psf with every exposure"
        );
    }
    // TODO: check if centroid is valid, unflagged
    afw::geom::Point2D position = measRecord.getCentroid();
    PTR(afw::detection::Psf::Image) psfImage = psf->computeImage(position);
    afw::geom::Box2I fitBBox = psfImage->getBBox();
    fitBBox.clip(exposure.getBBox());
    if (fitBBox != psfImage->getBBox()) {
        _flagHandler.setValue(measRecord, FAILURE, true);  // if we had a suspect flag, we'd set that instead
        _flagHandler.setValue(measRecord, EDGE, true);
    }
    afw::detection::Footprint fitRegion(fitBBox);
    if (!_ctrl.badMaskPlanes.empty()) {
        afw::image::MaskPixel badBits = 0x0;
        for (
            std::vector<std::string>::const_iterator i = _ctrl.badMaskPlanes.begin();
            i != _ctrl.badMaskPlanes.end();
            ++i
        ) {
            badBits |= exposure.getMaskedImage().getMask()->getPlaneBitMask(*i);
        }
        fitRegion.intersectMask(*exposure.getMaskedImage().getMask(), badBits);
    }
    if (fitRegion.getArea() == 0) {
        throw LSST_EXCEPT(
            MeasurementError,
            _flagHandler.getDefinition(NO_GOOD_PIXELS).doc,
            NO_GOOD_PIXELS
        );
    }
    typedef afw::detection::Psf::Pixel PsfPixel;
    typedef typename afw::image::MaskedImage<float>::Variance::Pixel VarPixel;
    ndarray::EigenView<PsfPixel,1,1,Eigen::ArrayXpr> model(
        afw::detection::flattenArray(
            fitRegion,
            psfImage->getArray(),
            psfImage->getXY0()
        )
    );
    ndarray::EigenView<float,1,1,Eigen::ArrayXpr> data(
        afw::detection::flattenArray(
            fitRegion,
            exposure.getMaskedImage().getImage()->getArray(),
            exposure.getXY0()
        )
    );
    ndarray::EigenView<VarPixel,1,1,Eigen::ArrayXpr> variance(
        afw::detection::flattenArray(
            fitRegion,
            exposure.getMaskedImage().getVariance()->getArray(),
            exposure.getXY0()
        )
    );
    PsfPixel alpha = model.matrix().squaredNorm();
    FluxResult result;
    result.flux = model.matrix().dot(data.matrix().cast<PsfPixel>()) / alpha;
    // If we're not using per-pixel weights to compute the flux, we'll still want to compute the
    // variance as if we had, so we'll apply the weights to the model vector now, and update alpha.
    result.fluxSigma = std::sqrt(model.square().matrix().dot(variance.matrix().cast<PsfPixel>()))
        / alpha;
    if (!utils::isfinite(result.flux) || !utils::isfinite(result.fluxSigma)) {
        throw LSST_EXCEPT(PixelValueError, "Invalid pixel value detected in image.");
    }
    measRecord.set(_fluxResultKey, result);
}

void PsfFluxAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}


}}} // namespace lsst::meas::base

