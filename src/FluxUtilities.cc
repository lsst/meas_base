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

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/Transform.h"

namespace lsst { namespace meas { namespace base {

FluxResult::FluxResult() :
    flux(std::numeric_limits<Flux>::quiet_NaN()),
    fluxSigma(std::numeric_limits<FluxErrElement>::quiet_NaN())
{}

FluxResultKey FluxResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc
) {
    FluxResultKey result;
    result._flux = schema.addField<Flux>(schema.join(name, "flux"), doc, "dn");
    result._fluxSigma = schema.addField<FluxErrElement>(schema.join(name, "fluxSigma"),
                                                        "1-sigma flux uncertainty", "dn");
    return result;
}

FluxResult FluxResultKey::get(afw::table::BaseRecord const & record) const {
    FluxResult r;
    r.flux = record.get(_flux);
    r.fluxSigma = record.get(_fluxSigma);
    return r;
}

void FluxResultKey::set(afw::table::BaseRecord & record, FluxResult const & value) const {
    record.set(_flux, value.flux);
    record.set(_fluxSigma, value.fluxSigma);
}

MagResultKey MagResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name
) {
    MagResultKey result;
    result._magKey = schema.addField<Mag>(schema.join(name, "mag"), "Magnitude");
    result._magErrKey = schema.addField<MagErrElement>(schema.join(name, "magErr"), "Error on magnitude");
    return result;
}

MagResult MagResultKey::get(afw::table::BaseRecord const & record) const {
    MagResult result = {record.get(_magKey), record.get(_magErrKey)};
    return result;
}

void MagResultKey::set(afw::table::BaseRecord & record, MagResult const & magResult) const {
    record.set(_magKey, magResult.mag);
    record.set(_magErrKey, magResult.magErr);
}

void MagResultKey::set(afw::table::BaseRecord & record, std::pair<double,double> const & magResult) const {
    record.set(_magKey, magResult.first);
    record.set(_magErrKey, magResult.second);
}

FluxTransform::FluxTransform(
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    BaseTransform{name}
{
    // Map the flag through to the output
    mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(name + "_flag").key);

    // Add keys for the magnitude and associated error
    _magKey = MagResultKey::addFields(mapper.editOutputSchema(), name);
}

void FluxTransform::operator()(
    afw::table::SourceCatalog const & inputCatalog,
    afw::table::BaseCatalog & outputCatalog,
    afw::image::Wcs const & wcs,
    afw::image::Calib const & calib
) const {
    checkCatalogSize(inputCatalog, outputCatalog);
    FluxResultKey fluxKey(inputCatalog.getSchema()[_name]);
    afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
    afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();
    {
        // While noThrow is in scope, converting a negative flux to a magnitude
        // returns NaN rather than throwing.
        NoThrowOnNegativeFluxContext noThrow;
        for (; inSrc < inputCatalog.end() && outSrc < outputCatalog.end(); ++inSrc, ++outSrc) {
            FluxResult fluxResult = fluxKey.get(*inSrc);
            _magKey.set(*outSrc, calib.getMagnitude(fluxResult.flux, fluxResult.fluxSigma));
        }
    }
}

NoThrowOnNegativeFluxContext::NoThrowOnNegativeFluxContext() {
    _throwOnNegative = afw::image::Calib::getThrowOnNegativeFlux();
    afw::image::Calib::setThrowOnNegativeFlux(false);
}

NoThrowOnNegativeFluxContext::~NoThrowOnNegativeFluxContext() {
    afw::image::Calib::setThrowOnNegativeFlux(_throwOnNegative);
}

}}} // namespace lsst::meas::base
