
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
#include "lsst/log/Log.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/meas/base/PsfFlux.h"

namespace lsst {
namespace meas {
namespace base {
namespace {
FlagDefinitionList flagDefinitions;
}  // namespace

FlagDefinition const PsfFluxAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
FlagDefinition const PsfFluxAlgorithm::NO_GOOD_PIXELS =
        flagDefinitions.add("flag_noGoodPixels", "not enough non-rejected pixels in data to attempt the fit");
FlagDefinition const PsfFluxAlgorithm::EDGE = flagDefinitions.add(
        "flag_edge", "object was too close to the edge of the image to use the full PSF model");

FlagDefinitionList const& PsfFluxAlgorithm::getFlagDefinitions() { return flagDefinitions; }

namespace {}  // namespace

PsfFluxAlgorithm::PsfFluxAlgorithm(Control const& ctrl, std::string const& name, afw::table::Schema& schema,
                                   std::string const& logName)
        : _ctrl(ctrl),
          _fluxResultKey(FluxResultKey::addFields(schema, name,
                                                  "flux derived from linear least-squares fit of PSF model")),
          _areaKey(schema.addField<float>(name + "_area", "effective area of PSF", "pixel")),
          _centroidExtractor(schema, name) {
    _logName = logName.size() ? logName : name;
    _flagHandler = FlagHandler::addFields(schema, name, getFlagDefinitions());
}

void PsfFluxAlgorithm::measure(afw::table::SourceRecord& measRecord,
                               afw::image::Exposure<float> const& exposure) const {
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        LOGL_ERROR(getLogName(), "PsfFlux: no psf attached to exposure");
        throw LSST_EXCEPT(FatalAlgorithmError, "PsfFlux algorithm requires a Psf with every exposure");
    }
    geom::Point2D position = _centroidExtractor(measRecord, _flagHandler);
    PTR(afw::detection::Psf::Image) psfImage = psf->computeImage(position);
    geom::Box2I fitBBox = psfImage->getBBox();
    fitBBox.clip(exposure.getBBox());
    if (fitBBox != psfImage->getBBox()) {
        _flagHandler.setValue(measRecord, FAILURE.number,
                              true);  // if we had a suspect flag, we'd set that instead
        _flagHandler.setValue(measRecord, EDGE.number, true);
    }
    auto fitRegionSpans = std::make_shared<afw::geom::SpanSet>(fitBBox);
    afw::detection::Footprint fitRegion(fitRegionSpans);
    if (!_ctrl.badMaskPlanes.empty()) {
        afw::image::MaskPixel badBits = 0x0;
        for (std::vector<std::string>::const_iterator i = _ctrl.badMaskPlanes.begin();
             i != _ctrl.badMaskPlanes.end(); ++i) {
            badBits |= exposure.getMaskedImage().getMask()->getPlaneBitMask(*i);
        }
        fitRegion.setSpans(fitRegion.getSpans()
                                   ->intersectNot(*exposure.getMaskedImage().getMask(), badBits)
                                   ->clippedTo(exposure.getMaskedImage().getMask()->getBBox()));
    }
    if (fitRegion.getArea() == 0) {
        throw LSST_EXCEPT(MeasurementError, NO_GOOD_PIXELS.doc, NO_GOOD_PIXELS.number);
    }
    typedef afw::detection::Psf::Pixel PsfPixel;
    auto model = fitRegion.getSpans()
                         ->flatten(psfImage->getArray(), psfImage->getXY0())
                         .asEigen<Eigen::ArrayXpr>();
    auto data = fitRegion.getSpans()
                        ->flatten(exposure.getMaskedImage().getImage()->getArray(), exposure.getXY0())
                        .asEigen<Eigen::ArrayXpr>();
    auto variance = fitRegion.getSpans()
                            ->flatten(exposure.getMaskedImage().getVariance()->getArray(), exposure.getXY0())
                            .asEigen<Eigen::ArrayXpr>();
    PsfPixel alpha = model.matrix().squaredNorm();
    FluxResult result;
    result.flux = model.matrix().dot(data.matrix().cast<PsfPixel>()) / alpha;
    // If we're not using per-pixel weights to compute the flux, we'll still want to compute the
    // variance as if we had, so we'll apply the weights to the model vector now, and update alpha.
    result.fluxSigma = std::sqrt(model.square().matrix().dot(variance.matrix().cast<PsfPixel>())) / alpha;
    measRecord.set(_areaKey, model.matrix().sum() / alpha);
    if (!std::isfinite(result.flux) || !std::isfinite(result.fluxSigma)) {
        throw LSST_EXCEPT(PixelValueError, "Invalid pixel value detected in image.");
    }
    measRecord.set(_fluxResultKey, result);
}

void PsfFluxAlgorithm::fail(afw::table::SourceRecord& measRecord, MeasurementError* error) const {
    _flagHandler.handleFailure(measRecord, error);
}

PsfFluxTransform::PsfFluxTransform(Control const& ctrl, std::string const& name,
                                   afw::table::SchemaMapper& mapper)
        : FluxTransform{name, mapper} {
    for (std::size_t i = 0; i < PsfFluxAlgorithm::getFlagDefinitions().size(); i++) {
        FlagDefinition const& flag = PsfFluxAlgorithm::getFlagDefinitions()[i];
        if (flag == PsfFluxAlgorithm::FAILURE) continue;
        if (mapper.getInputSchema().getNames().count(mapper.getInputSchema().join(name, flag.name)) == 0)
            continue;
        afw::table::Key<afw::table::Flag> key =
                mapper.getInputSchema().find<afw::table::Flag>(name + "_" + flag.name).key;
        mapper.addMapping(key);
    }
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
