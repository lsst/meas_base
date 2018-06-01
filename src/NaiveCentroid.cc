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
#include "ndarray/eigen.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/NaiveCentroid.h"

namespace lsst {
namespace meas {
namespace base {
namespace {
FlagDefinitionList flagDefinitions;
}  // namespace

FlagDefinition const NaiveCentroidAlgorithm::FAILURE = flagDefinitions.addFailureFlag();
FlagDefinition const NaiveCentroidAlgorithm::NO_COUNTS =
        flagDefinitions.add("flag_noCounts", "Object to be centroided has no counts");
FlagDefinition const NaiveCentroidAlgorithm::EDGE =
        flagDefinitions.add("flag_edge", "Object too close to edge");

FlagDefinitionList const& NaiveCentroidAlgorithm::getFlagDefinitions() { return flagDefinitions; }

NaiveCentroidAlgorithm::NaiveCentroidAlgorithm(Control const& ctrl, std::string const& name,
                                               afw::table::Schema& schema)
        : _ctrl(ctrl),
          _centroidKey(CentroidResultKey::addFields(schema, name, "centroid from Naive Centroid algorithm",
                                                    NO_UNCERTAINTY)),
          _flagHandler(FlagHandler::addFields(schema, name, getFlagDefinitions())),
          _centroidExtractor(schema, name, true),
          _centroidChecker(schema, name, ctrl.doFootprintCheck, ctrl.maxDistToPeak) {}

void NaiveCentroidAlgorithm::measure(afw::table::SourceRecord& measRecord,
                                     afw::image::Exposure<float> const& exposure) const {
    geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    CentroidResult result;
    result.x = center.getX();
    result.y = center.getY();
    measRecord.set(_centroidKey, result);  // better than NaN

    typedef afw::image::Image<float> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = center.getX();  // FIXME: this is different from GaussianCentroid and SdssCentroid here,
    int y = center.getY();  //        and probably shouldn't be.

    x -= image.getX0();  // work in image Pixel coordinates
    y -= image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {
        throw LSST_EXCEPT(MeasurementError, EDGE.doc, EDGE.number);
    }

    ImageT::xy_locator im = image.xy_at(x, y);

    double const sum = (im(-1, 1) + im(0, 1) + im(1, 1) + im(-1, 0) + im(0, 0) + im(1, 0) + im(-1, -1) +
                        im(0, -1) + im(1, -1)) -
                       9 * _ctrl.background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(MeasurementError, NO_COUNTS.doc, NO_COUNTS.number);
    }

    double const sum_x = -im(-1, 1) + im(1, 1) + -im(-1, 0) + im(1, 0) + -im(-1, -1) + im(1, -1);
    double const sum_y = (im(-1, 1) + im(0, 1) + im(1, 1)) - (im(-1, -1) + im(0, -1) + im(1, -1));

    result.x = afw::image::indexToPosition(x + image.getX0()) + sum_x / sum;
    result.y = afw::image::indexToPosition(y + image.getY0()) + sum_y / sum;
    measRecord.set(_centroidKey, result);
    _centroidChecker(measRecord);
}

void NaiveCentroidAlgorithm::fail(afw::table::SourceRecord& measRecord, MeasurementError* error) const {
    _flagHandler.handleFailure(measRecord, error);
}

NaiveCentroidTransform::NaiveCentroidTransform(Control const& ctrl, std::string const& name,
                                               afw::table::SchemaMapper& mapper)
        : CentroidTransform{name, mapper} {
    for (std::size_t i = 0; i < NaiveCentroidAlgorithm::getFlagDefinitions().size(); i++) {
        FlagDefinition const& flag = NaiveCentroidAlgorithm::getFlagDefinitions()[i];
        if (flag == NaiveCentroidAlgorithm::FAILURE) continue;
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
