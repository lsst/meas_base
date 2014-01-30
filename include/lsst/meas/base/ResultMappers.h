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
template <std::size_t N>
struct FlagsResultMapper {

    FlagsResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        boost::array<FlagDef,N> const & flagDefs
    );

    // Transfer values from the result struct to the record, and set the failure flag field.
    void fail(afw::table::BaseRecord & record) const;

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, FlagsResult<N> const & result) const;

protected:
    boost::array<afw::table::Key<afw::table::Flag>,N+1> _flags;
};

// Object that transfers values from FluxResultResult to a record object.
class FluxResultMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    FluxResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, FluxResult const & result) const;

private:
    afw::table::Key<Flux> _flux;
    afw::table::Key<ErrElement> _fluxSigma;
};

// Object that transfers values from CentroidResultResult to a record object.
class CentroidResultMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    CentroidResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, CentroidResult const & result) const;

private:
    afw::table::Key<CentroidElement> _x;
    afw::table::Key<CentroidElement> _y;
    afw::table::Key<ErrElement> _xSigma;
    afw::table::Key<ErrElement> _ySigma;
    afw::table::Key<ErrElement> _x_y_Cov;
};

// Object that transfers values from ShapeResultResult to a record object.
class ShapeResultMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    ShapeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, ShapeResult const & result) const;

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

template <typename Algorithm, typename T1>
struct SimpleResultMapper1 : public T1, public FlagsResultMapper<Algorithm::N_FLAGS> {

    SimpleResultMapper1(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1
    ) : T1(schema, prefix, uncertainty1),
        FlagsResultMapper<Algorithm::N_FLAGS>(schema, prefix, Algorithm::getFlagDefinitions())
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        FlagsResultMapper<Algorithm::N_FLAGS>(record, result);
    }

};

template <typename Algorithm, typename T1, typename T2>
struct SimpleResultMapper2 : public T1, public T2, public FlagsResultMapper<Algorithm::N_FLAGS> {

    SimpleResultMapper2(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1,
        ResultMapperUncertaintyEnum uncertainty2
    ) : T1(schema, prefix, uncertainty1),
        T2(schema, prefix, uncertainty2),
        FlagsResultMapper<Algorithm::N_FLAGS>(schema, prefix, Algorithm::getFlagDefinitions())
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        T2::apply(record, result);
        FlagsResultMapper<Algorithm::N_FLAGS>(record, result);
    }

};

template <typename Algorithm, typename T1, typename T2, typename T3>
struct SimpleResultMapper3 : public T1, public T2, public T3, public FlagsResultMapper<Algorithm::N_FLAGS> {

    SimpleResultMapper3(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1,
        ResultMapperUncertaintyEnum uncertainty2,
        ResultMapperUncertaintyEnum uncertainty3
    ) : T1(schema, prefix, uncertainty1),
        T2(schema, prefix, uncertainty2),
        T3(schema, prefix, uncertainty3),
        FlagsResultMapper<Algorithm::N_FLAGS>(schema, prefix, Algorithm::getFlagDefinitions())
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        T2::apply(record, result);
        T3::apply(record, result);
        FlagsResultMapper<Algorithm::N_FLAGS>(record, result);
    }

};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ResultMappers_h_INCLUDED
