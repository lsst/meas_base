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
#include "lsst/utils/ieee.h"

namespace lsst { namespace meas { namespace base {

FootprintInput::Vector FootprintInput::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(FootprintInput(*i));
    }
    return r;
}

FootprintCentroidInput::FootprintCentroidInput(afw::table::SourceRecord const & record) : FootprintInput(record) 
{
        
        if (!record.getTable()->getCentroidKey().isValid() || lsst::utils::isnan(record.getCentroid().getX())) {
            position.setX(record.getFootprint()->getPeaks()[0]->getFx());
            position.setY(record.getFootprint()->getPeaks()[0]->getFy());
        } else {
            position = record.getCentroid();
        }
}

FootprintCentroidInput::Vector FootprintCentroidInput::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(FootprintCentroidInput(*i));
    }
    return r;
}

FootprintCentroidShapeInput::Vector FootprintCentroidShapeInput::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(FootprintCentroidShapeInput(*i));
    }
    return r;
}

}}} // namespace lsst::meas::base
