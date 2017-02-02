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

#include <array>
#include <cmath>

#include "ndarray/eigen.h"

#include "lsst/afw/table/Source.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/meas/base/PsfFlux.h"

namespace lsst { namespace meas { namespace base {
struct PsfFluxAlgorithm::Flags {
    static FlagDefinition FAILURE;
    static FlagDefinition NO_GOOD_PIXELS;
    static FlagDefinition EDGE;
};
FlagDefinition PsfFluxAlgorithm::Flags::FAILURE("flag", "general failure flag");
FlagDefinition PsfFluxAlgorithm::Flags::NO_GOOD_PIXELS("flag_noGoodPixels", "not enough non-rejected pixels in data to attempt the fit");
FlagDefinition PsfFluxAlgorithm::Flags::EDGE("flag_edge", "object was too close to the edge of the image to use the full PSF model");
namespace {
std::vector<FlagDefinition> const flagVector = {
    PsfFluxAlgorithm::Flags::FAILURE,
    PsfFluxAlgorithm::Flags::NO_GOOD_PIXELS,
    PsfFluxAlgorithm::Flags::EDGE,
};
std::vector<FlagDefinition> const & getFlagDefinitions() {
    return flagVector;
};
} // end anonymous

std::size_t PsfFluxAlgorithm::getFlagNumber(std::string const & name) {
    std::size_t i = 0;
    for (auto iter = getFlagDefinitions().begin(); iter < getFlagDefinitions().end(); iter++) {
        if (iter->name == name) {
            return i;
        }
        i++;
    }
    throw lsst::pex::exceptions::RuntimeError("PsfFlux flag does not exist for name: " + name);
}

std::string const PsfFluxAlgorithm::getFlagName(std::size_t flagNumber) {
    std::size_t i = 0;
    for (auto iter = getFlagDefinitions().begin(); iter < getFlagDefinitions().end(); iter++) {
        if (i == flagNumber) {
            return iter->name;
        }
    }
    throw lsst::pex::exceptions::RuntimeError("PsfFlux flag does not exist for number: " + flagNumber);
}

namespace {


} // anonymous

PsfFluxAlgorithm::PsfFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxResultKey(
        FluxResultKey::addFields(schema, name, "flux derived from linear least-squares fit of PSF model")
    ),
    _centroidExtractor(schema, name)
{
    _flagHandler = FlagHandler::addFields(schema, name,
                                          getFlagDefinitions().begin(), getFlagDefinitions().end());
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
    afw::geom::Point2D position = _centroidExtractor(measRecord, _flagHandler);
    PTR(afw::detection::Psf::Image) psfImage = psf->computeImage(position);
    afw::geom::Box2I fitBBox = psfImage->getBBox();
    fitBBox.clip(exposure.getBBox());
    if (fitBBox != psfImage->getBBox()) {
        _flagHandler.setValue(measRecord, Flags::FAILURE.name, true);  // if we had a suspect flag, we'd set that instead
        _flagHandler.setValue(measRecord, Flags::EDGE.name, true);
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
            Flags::NO_GOOD_PIXELS.doc,
            getFlagNumber(Flags::NO_GOOD_PIXELS.name)
        );
    }
    typedef afw::detection::Psf::Pixel PsfPixel;
    typedef afw::image::MaskedImage<float>::Variance::Pixel VarPixel;
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
    if (!std::isfinite(result.flux) || !std::isfinite(result.fluxSigma)) {
        throw LSST_EXCEPT(PixelValueError, "Invalid pixel value detected in image.");
    }
    measRecord.set(_fluxResultKey, result);
}

void PsfFluxAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

PsfFluxTransform::PsfFluxTransform(
    Control const & ctrl,
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    FluxTransform{name, mapper}
{
    for (auto flag = getFlagDefinitions().begin() + 1; flag < getFlagDefinitions().end(); flag++) {
        mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(name + "_" + flag->name).key);
    }
}

}}} // namespace lsst::meas::base
