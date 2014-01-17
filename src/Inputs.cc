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

#include "lsst/meas/base/Inputs.h"

namespace lsst { namespace meas { namespace base {

AlgorithmInput1::Vector AlgorithmInput1::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(AlgorithmInput1(*i));
    }
    return r;
}

AlgorithmInput2::Vector AlgorithmInput2::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(AlgorithmInput2(*i));
    }
    return r;
}

AlgorithmInput3::Vector AlgorithmInput3::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(AlgorithmInput3(*i));
    }
    return r;
}

}}} // namespace lsst::meas::base
