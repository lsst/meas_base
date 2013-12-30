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

#ifndef LSST_MEAS_BASE_ResultMappers_h_INCLUDED
#define LSST_MEAS_BASE_ResultMappers_h_INCLUDED

#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/Results.h"

namespace lsst { namespace meas { namespace base {

// Base class for classes that transfer values from Result structs to record objects
struct BaseAlgorithmMapper {

    BaseAlgorithmMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    // Set the failure flag field.
    void fail(afw::table::BaseRecord & record);

protected:
    afw::table::Key<afw::table::Flag> _flag;
};

// Object that transfers values from FluxAlgorithmResult to a record object.
class FluxAlgorithmMapper : public BaseAlgorithmMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    FluxAlgorithmMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, FluxAlgorithmResult const & result);

private:
    afw::table::Key<Flux> _value;
    afw::table::Key<FluxErr> _err;
};

// Object that transfers values from CentroidAlgorithmResult to a record object.
class CentroidAlgorithmMapper : public BaseAlgorithmMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    CentroidAlgorithmMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, CentroidAlgorithmResult const & result);

private:
    // We should be able to derive these Key types from the typedefs at the top of the file,
    // but that's a problem for afw::table to solve.
    afw::table::Key< afw::table::Point<double> > _value;
    afw::table::Key< afw::table::Covariance< afw::table::Point<float> > > _cov;
};

// Object that transfers values from ShapeAlgorithmResult to a record object.
class ShapeAlgorithmMapper : public BaseAlgorithmMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    ShapeAlgorithmMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, ShapeAlgorithmResult const & result);

private:
    // We should be able to derive these Key types from the typedefs at the top of the file,
    // but that's a problem for afw::table to solve.
    afw::table::Key< afw::table::Moments<double> > _value;
    afw::table::Key< afw::table::Covariance< afw::table::Moments<float> > > _cov;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ResultMappers_h_INCLUDED
