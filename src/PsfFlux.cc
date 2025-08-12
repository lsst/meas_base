
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
          _instFluxResultKey(FluxResultKey::addFields(
                  schema, name, "instFlux derived from linear least-squares fit of PSF model")),
          _areaKey(schema.addField<float>(name + "_area", "effective area of PSF", "pixel")),
          _chi2Key(schema.addField<float>(name + "_chi2", "chi2 of the fitted PSF")),
          _npixelsKey(schema.addField<int>(name + "_npixels",
                                           "number of pixels that were included in the PSF fit", "pixel")),
          _centroidExtractor(schema, name) {
    _logName = logName.size() ? logName : name;
    _flagHandler = FlagHandler::addFields(schema, name, getFlagDefinitions());
}

void PsfFluxAlgorithm::measure(afw::table::SourceRecord& measRecord,
                               afw::image::Exposure<float> const& exposure) const {
    std::shared_ptr<afw::detection::Psf const> psf = exposure.getPsf();
    if (!psf) {
        LOGL_ERROR(getLogName(), "PsfFlux: no psf attached to exposure");
        throw LSST_EXCEPT(FatalAlgorithmError, "PsfFlux algorithm requires a Psf with every exposure");
    }
    geom::Point2D position = _centroidExtractor(measRecord, _flagHandler);
    std::shared_ptr<afw::detection::Psf::Image> psfImage = psf->computeImage(position);
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
        _flagHandler.setValue(measRecord, NO_GOOD_PIXELS.number, true);
        _flagHandler.setValue(measRecord, FAILURE.number, true);
        return;
    }
    typedef afw::detection::Psf::Pixel PsfPixel;
    // SpanSet::flatten returns a new ndarray::Array, which must stay in scope
    // while we use an Eigen::Map view of it
    auto modelNdArray = fitRegion.getSpans()->flatten(psfImage->getArray(), psfImage->getXY0());
    auto dataNdArray = fitRegion.getSpans()->flatten(exposure.getMaskedImage().getImage()->getArray(),
                                                     exposure.getXY0());
    auto varianceNdArray = fitRegion.getSpans()->flatten(exposure.getMaskedImage().getVariance()->getArray(),
                                                         exposure.getXY0());
    auto model = ndarray::asEigenMatrix(modelNdArray);
    auto data = ndarray::asEigenMatrix(dataNdArray);
    auto variance = ndarray::asEigenMatrix(varianceNdArray);
    PsfPixel alpha = model.squaredNorm();
    FluxResult result;
    result.instFlux = model.dot(data.cast<PsfPixel>()) / alpha;
    // If we're not using per-pixel weights to compute the instFlux, we'll still want to compute the
    // variance as if we had, so we'll apply the weights to the model now, and update alpha.
    result.instFluxErr = std::sqrt(model.array().square().matrix().dot(variance.cast<PsfPixel>())) / alpha;
    measRecord.set(_areaKey, std::pow(model.sum(), 2) / alpha);
    measRecord.set(_npixelsKey, fitRegion.getSpans()->getArea());
    if (!std::isfinite(result.instFlux) || !std::isfinite(result.instFluxErr)) {
        throw LSST_EXCEPT(PixelValueError, "Invalid pixel value detected in image.");
    }
    measRecord.set(_instFluxResultKey, result);
    auto chi2 = ((data.cast<PsfPixel>() - result.instFlux * model).array().square() /
                 variance.cast<PsfPixel>().array())
                        .sum();
    measRecord.set(_chi2Key, chi2);
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
