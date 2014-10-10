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

#ifndef LSST_MEAS_BASE_FlagHandler_h_INCLUDED
#define LSST_MEAS_BASE_FlagHandler_h_INCLUDED

#include <vector>

#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/exceptions.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief Simple POD struct used to define and document flags
 *
 *  When we switch to C++11, we can make the attributes std::strings, but at least for now we'll use
    C-strings so we can create these arrays using initializer lists even in C++98.
 */
struct FlagDefinition {
    char const * name;
    char const * doc;
};

class FlagHandler {
public:

    FlagHandler() {}

    explicit FlagHandler(
        afw::table::Schema & schema,
        std::string const & prefix,
        FlagDefinition const * begin,
        FlagDefinition const * end
    );

    FlagDefinition getDefinition(int i) const { return _vector[i].first; }

    bool getValue(afw::table::BaseRecord const & record, int i) const {
        return record.get(_vector[i].second);
    }

    void setValue(afw::table::BaseRecord & record, int i, bool value) const {
        record.set(_vector[i].second, value);
    }

    void handleFailure(afw::table::BaseRecord & record, MeasurementError const * error=NULL) const;

private:

    typedef std::vector< std::pair<FlagDefinition, afw::table::Key<afw::table::Flag> > > Vector;

    Vector _vector;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_FlagHandler_h_INCLUDED
