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

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/ScaledApertureFlux.h"
#include "lsst/meas/base/SincCoeffs.h"
#include "lsst/afw/detection/Psf.h"

namespace lsst { namespace meas { namespace base {

ScaledApertureFluxAlgorithm::ScaledApertureFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxResultKey(
        FluxResultKey::addFields(schema, name, "flux derived from PSF-scaled aperture")
    ),
    _centroidExtractor(schema, name)
{
    _flagHandler = FlagHandler::addFields(schema, name,
                                          ApertureFluxAlgorithm::getFlagDefinitions().begin(),
                                          ApertureFluxAlgorithm::getFlagDefinitions().end());
}

void ScaledApertureFluxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    afw::geom::Point2D const center = _centroidExtractor(measRecord, _flagHandler);
    double const radius = exposure.getPsf()->computeShape(center).getDeterminantRadius();
    double const fwhm = 2.0*std::sqrt(2.0*std::log(2))*radius;
    double const size = _ctrl.scale*fwhm;
    afw::geom::ellipses::Axes const axes(size, size);

    // ApertureFluxAlgorithm::computeSincFlux requires an ApertureFluxControl as an
    // argument. All that it uses it for is to read the type of warping kernel.
    ApertureFluxControl apCtrl;
    apCtrl.shiftKernel = _ctrl.shiftKernel;

    Result result = ApertureFluxAlgorithm::computeSincFlux(exposure.getMaskedImage(),
                                                           afw::geom::ellipses::Ellipse(axes, center),
                                                           apCtrl);
    measRecord.set(_fluxResultKey, result);

    for (FlagDefinition const * iter = ApertureFluxAlgorithm::getFlagDefinitions().begin();
         iter < ApertureFluxAlgorithm::getFlagDefinitions().end(); iter++) {
        if (result.getFlag(iter->number)) {
            _flagHandler.setValue(measRecord, iter->number, true);
        }
    }
}

void ScaledApertureFluxAlgorithm::fail(afw::table::SourceRecord & measRecord,
                                       MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

ScaledApertureFluxTransform::ScaledApertureFluxTransform(
    Control const & ctrl,
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    FluxTransform{name, mapper}
{
    for (auto flag = ApertureFluxAlgorithm::getFlagDefinitions().begin() + 1;
         flag < ApertureFluxAlgorithm::getFlagDefinitions().end(); flag++) {
        mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(name + "_" + flag->name).key);
    }
}

}}} // namespace lsst::meas::base
