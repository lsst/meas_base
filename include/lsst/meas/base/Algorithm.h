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

#ifndef LSST_MEAS_BASE_Algorithm_h_INCLUDED
#define LSST_MEAS_BASE_Algorithm_h_INCLUDED

#include "lsst/afw/table/fwd.h"
#include "lsst/meas/base/exceptions.h"

namespace lsst { namespace meas { namespace base {

class BaseAlgorithm {
public:

    void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const = 0;

    virtual ~BaseAlgorithm() {}

};

class SingleFrameAlgorithm : public virtual BaseAlgorithm {
public:

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const = 0;

    void measureN(
        afw::table::SourceCatalog const & measCat,
        afw::image::Exposure<float> const & exposure
    ) const;

};

class ForcedAlgorithm : public virtual BaseAlgorithm {
public:

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceRecord const & refRecord,
        afw::image::Wcs const & refWcs
    ) const = 0;

    void measureN(
        afw::table::SourceCatalog const & measCat,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceCatalog const & refRecord,
        afw::image::Wcs const & refWcs
    ) const;

};

class SimpleAlgorithm : public SingleFrameAlgorithm, public ForcedAlgorithm {
public:

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceRecord const & refRecord,
        afw::image::Wcs const & refWcs
    ) const {
        measure(measRecord, exposure);
    }

    void measureN(
        afw::table::SourceCatalog const & measCat,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceCatalog const & refRecord,
        afw::image::Wcs const & refWcs
    ) const {
        measureN(measCat, exposure);
    }

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_Algorithm_h_INCLUDED
