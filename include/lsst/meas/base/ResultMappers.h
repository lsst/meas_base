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

enum ResultMapperUncertaintyEnum {
    NO_UNCERTAINTY = 0,
    DIAGONAL_ONLY = 1,
    FULL_COVARIANCE = 2
};

// Base class for classes that transfer values from Result structs to record objects
struct BaseResultMapper {

    BaseResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    // Set the failure flag field.
    void fail(afw::table::BaseRecord & record);

protected:
    afw::table::Key<afw::table::Flag> _flag;
};

// Object that transfers values from FluxResultResult to a record object.
class FluxResultMapper : public BaseResultMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    FluxResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, FluxResult const & result);

private:
    afw::table::Key<Flux> _flux;
    afw::table::Key<ErrElement> _fluxSigma;
};

// Object that transfers values from CentroidResultResult to a record object.
class CentroidResultMapper : public BaseResultMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    CentroidResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, CentroidResult const & result);

private:
    afw::table::Key<CentroidElement> _x;
    afw::table::Key<CentroidElement> _y;
    afw::table::Key<ErrElement> _xSigma;
    afw::table::Key<ErrElement> _ySigma;
    afw::table::Key<ErrElement> _x_y_Cov;
};

// Object that transfers values from ShapeResultResult to a record object.
class ShapeResultMapper : public BaseResultMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    ShapeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, ShapeResult const & result);

private:
    afw::table::Key<ShapeElement> _xx;
    afw::table::Key<ShapeElement> _yy;
    afw::table::Key<ShapeElement> _xy;
    afw::table::Key<ErrElement> _xxSigma;
    afw::table::Key<ErrElement> _yySigma;
    afw::table::Key<ErrElement> _xySigma;
    afw::table::Key<ErrElement> _xx_yy_Cov;
    afw::table::Key<ErrElement> _xx_xy_Cov;
    afw::table::Key<ErrElement> _yy_xy_Cov;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ResultMappers_h_INCLUDED
