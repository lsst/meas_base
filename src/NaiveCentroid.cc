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


namespace lsst { namespace meas { namespace base {

namespace {

std::array<FlagDefinition,NaiveCentroidAlgorithm::N_FLAGS> const & getFlagDefinitions() {
    static std::array<FlagDefinition,NaiveCentroidAlgorithm::N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag, set if anything went wrong"},
        {"flag_noCounts", "Object to be centroided has no counts"},
        {"flag_edge", "Object too close to edge"}
    }};
    return flagDefs;
}

} // anonymous

NaiveCentroidAlgorithm::NaiveCentroidAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _centroidKey(
        CentroidResultKey::addFields(schema, name, "centroid from Naive Centroid algorithm", NO_UNCERTAINTY)
    ),
    _centroidExtractor(schema, name, true)
{
    _flagHandler = FlagHandler::addFields(schema, name,
                                          getFlagDefinitions().begin(), getFlagDefinitions().end());
}

void NaiveCentroidAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    
    afw::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    CentroidResult result;
    result.x = center.getX();
    result.y = center.getY(); 
    measRecord.set(_centroidKey, result); // better than NaN

    typedef afw::image::Image<float> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = center.getX();  // FIXME: this is different from GaussianCentroid and SdssCentroid here,
    int y = center.getY();  //        and probably shouldn't be.

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {

        throw LSST_EXCEPT(
            MeasurementError,
            _flagHandler.getDefinition(EDGE).doc,
            EDGE
        );
    }

    ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1))
        - 9 * _ctrl.background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(
            MeasurementError,
            _flagHandler.getDefinition(NO_COUNTS).doc,
            NO_COUNTS
        );
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    result.x = lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x / sum;
    result.y = lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y / sum;
    measRecord.set(_centroidKey, result);
    _flagHandler.setValue(measRecord, FAILURE, false);  // if we had a suspect flag, we'd set that instead
}


void NaiveCentroidAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

NaiveCentroidTransform::NaiveCentroidTransform(
    Control const & ctrl,
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    CentroidTransform{name, mapper}
{
    for (auto flag = getFlagDefinitions().begin() + 1; flag < getFlagDefinitions().end(); ++flag) {
        mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(
                          mapper.getInputSchema().join(name, flag->name)).key);
    }
}

}}} // namespace lsst::meas::base


