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

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/GaussianFlux.h"
#include "lsst/meas/base/SdssShape.h"

namespace lsst { namespace meas { namespace base {
struct GaussianFluxAlgorithm::Flags {
    static FlagDefinition FAILURE;
};
FlagDefinition GaussianFluxAlgorithm::Flags::FAILURE("flag", "general failure flag, set if anything went wrong");
namespace {
std::vector<FlagDefinition> const flagVector = {
    GaussianFluxAlgorithm::Flags::FAILURE,
};
std::vector<FlagDefinition> const & getFlagDefinitions() {
    return flagVector;
};
} // end anonymous

std::size_t GaussianFluxAlgorithm::getFlagNumber(std::string const & name) {
    std::size_t i = 0;
    for (auto iter = getFlagDefinitions().begin(); iter < getFlagDefinitions().end(); iter++) {
        if (iter->name == name) {
            return i;
        }
        i++;
    }
    throw lsst::pex::exceptions::RuntimeError("GaussianFlux flag does not exist for name: " + name);
}

std::string const GaussianFluxAlgorithm::getFlagName(std::size_t flagNumber) {
    std::size_t i = 0;
    for (auto iter = getFlagDefinitions().begin(); iter < getFlagDefinitions().end(); iter++) {
        if (i == flagNumber) {
            return iter->name;
        }
    }
    throw lsst::pex::exceptions::RuntimeError("GaussianFlux flag does not exist for number: " + flagNumber);
}

GaussianFluxAlgorithm::GaussianFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxResultKey(
        FluxResultKey::addFields(schema, name, "flux from Gaussian Flux algorithm")
    ),
    _centroidExtractor(schema, name),
    _shapeExtractor(schema, name)
{
    _flagHandler = FlagHandler::addFields(schema, name, getFlagDefinitions().begin(), getFlagDefinitions().end());
}

void GaussianFluxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    afw::geom::Point2D centroid = _centroidExtractor(measRecord, _flagHandler);
    afw::geom::ellipses::Quadrupole shape = _shapeExtractor(measRecord, _flagHandler);

    FluxResult result = SdssShapeAlgorithm::computeFixedMomentsFlux(
        exposure.getMaskedImage(), shape, centroid
    );

    measRecord.set(_fluxResultKey, result);
    _flagHandler.setValue(measRecord, Flags::FAILURE.name, false);
}


void GaussianFluxAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}} // namespace lsst::meas::base
