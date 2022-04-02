// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018 AURA/LSST.
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

#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/table/Source.h"
#include "lsst/log/Log.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/base/LocalBackground.h"

namespace lsst {
namespace meas {
namespace base {
namespace {
FlagDefinitionList flagDefinitions;
}  // namespace

FlagDefinition const LocalBackgroundAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
FlagDefinition const LocalBackgroundAlgorithm::NO_GOOD_PIXELS =
        flagDefinitions.add("flag_noGoodPixels", "no good pixels in the annulus");
FlagDefinition const LocalBackgroundAlgorithm::NO_PSF = flagDefinitions.add("flag_noPsf", "no PSF provided");

FlagDefinitionList const& LocalBackgroundAlgorithm::getFlagDefinitions() { return flagDefinitions; }

LocalBackgroundAlgorithm::LocalBackgroundAlgorithm(Control const& ctrl, std::string const& name,
                                                   afw::table::Schema& schema, std::string const& logName)
        : _ctrl(ctrl),
          _resultKey(FluxResultKey::addFields(schema, name, "background in annulus around source")),
          _flagHandler(FlagHandler::addFields(schema, name, getFlagDefinitions())),
          _centroidExtractor(schema, name),
          _stats(ctrl.bgRej, ctrl.bgIter) {
    _logName = logName.size() ? logName : name;
}

void LocalBackgroundAlgorithm::measure(afw::table::SourceRecord& measRecord,
                                       afw::image::Exposure<float> const& exposure) const {
    geom::Point2D const center = _centroidExtractor(measRecord, _flagHandler);
    afw::image::MaskedImage<float> const& image = exposure.getMaskedImage();
    afw::image::Mask<afw::image::MaskPixel> const& mask = *image.getMask();
    afw::image::MaskPixel const badMask = mask.getPlaneBitMask(_ctrl.badMaskPlanes);

    // Define pixels in annulus
    auto const psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(MeasurementError, NO_PSF.doc, NO_PSF.number);
    }
    float const psfSigma = psf->computeShape().getDeterminantRadius();

    float const innerRadius = _ctrl.annulusInner * psfSigma;
    afw::geom::ellipses::Axes const innerCircle{innerRadius, innerRadius};
    auto const& inner = afw::geom::SpanSet::fromShape(afw::geom::ellipses::Ellipse(innerCircle, center));

    float const outerRadius = _ctrl.annulusOuter * psfSigma;
    afw::geom::ellipses::Axes const outerCircle{outerRadius, outerRadius};
    auto const& outer = afw::geom::SpanSet::fromShape(afw::geom::ellipses::Ellipse(outerCircle, center));

    auto const& annulus = outer->clippedTo(image.getBBox())->intersectNot(*inner);
    auto const& imageValues = annulus->flatten(image.getImage()->getArray(), image.getXY0());
    auto const& maskValues = annulus->flatten(image.getMask()->getArray(), image.getXY0());

    // Extract from ndarray::Array into std::vector because of limitations in afw::math::makeStatistics
    std::vector<float> values;
    values.reserve(imageValues.getNumElements());
    if(imageValues.getNumElements() != maskValues.getNumElements()) {
        throw LSST_EXCEPT(MeasurementError, NO_GOOD_PIXELS.doc, NO_GOOD_PIXELS.number);
    }  // constructed from the same spans
    auto maskIter = maskValues.begin();
    for (auto imageIter = imageValues.begin(); imageIter != imageValues.end(); ++imageIter, ++maskIter) {
        if ((*maskIter & badMask) == 0) {
            values.push_back(*imageIter);
        }
    }

    if (values.size() == 0) {
        throw LSST_EXCEPT(MeasurementError, NO_GOOD_PIXELS.doc, NO_GOOD_PIXELS.number);
    }

    // Measure the background
    auto const stats = afw::math::makeStatistics(values, afw::math::MEANCLIP | afw::math::STDEVCLIP, _stats);
    FluxResult const result(stats.getValue(afw::math::MEANCLIP), stats.getValue(afw::math::STDEVCLIP));
    measRecord.set(_resultKey, result);
}

void LocalBackgroundAlgorithm::fail(afw::table::SourceRecord& measRecord, MeasurementError* error) const {
    _flagHandler.handleFailure(measRecord, error);
}

LocalBackgroundTransform::LocalBackgroundTransform(Control const& ctrl, std::string const& name,
                                                   afw::table::SchemaMapper& mapper)
        : FluxTransform{name, mapper} {
    for (std::size_t i = 0; i < LocalBackgroundAlgorithm::getFlagDefinitions().size(); i++) {
        FlagDefinition const& flag = LocalBackgroundAlgorithm::getFlagDefinitions()[i];
        if (flag == LocalBackgroundAlgorithm::FAILURE) continue;
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
