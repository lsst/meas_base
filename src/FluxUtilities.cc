// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/afw/table/BaseRecord.h"

namespace lsst { namespace meas { namespace base {

FluxResult::FluxResult() :
    flux(std::numeric_limits<Flux>::quiet_NaN()),
    fluxSigma(std::numeric_limits<FluxErrElement>::quiet_NaN())
{}

FluxResultKey FluxResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc,
    UncertaintyEnum uncertainty

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

}}} // namespace lsst::meas::base
