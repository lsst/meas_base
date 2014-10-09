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

#include "lsst/meas/base/FlagHandler.h"

namespace lsst { namespace meas { namespace base {

FlagHandler::FlagHandler(
    afw::table::Schema & schema,
    std::string const & prefix,
    FlagDefinition const * begin,
    FlagDefinition const * end
) {
    _keys.reserve(end - begin);
    for (FlagDefinition const * iter = begin; iter != end; ++iter) {
        _keys.push_back(
            schema.addField<afw::table::Flag>(schema.join(prefix, iter->name), iter->doc)
        );
    }
}

void FlagHandler::fail(afw::table::BaseRecord & record, MeasurementError const * error) const {
    record.set(_keys[0], true);
    if (error) {
        record.set(_keys[error->getFlagBit()], true);
    }
}

}}} // lsst::meas::base
