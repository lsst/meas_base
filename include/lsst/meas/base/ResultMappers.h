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

template <typename Algorithm>
struct FlagsComponentMapper {

    FlagsComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    // Set the general failure flag field and the specific bit held by the given exception.
    void fail(afw::table::BaseRecord & record, MeasurementError const * error=NULL) const;

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, FlagsComponent<Algorithm> const & result) const;

protected:
    boost::array<afw::table::Key<afw::table::Flag>,Algorithm::N_FLAGS+1> _flags;
};

template <typename Algorithm>
FlagsComponentMapper<Algorithm>::FlagsComponentMapper(
    afw::table::Schema & schema,
    std::string const & prefix
) {
    _flags[0] = schema.addField(
        afw::table::Field<afw::table::Flag>(
            prefix + "_flag",
            "general failure flag for " + prefix + " measurement"
        ),
        true // replace existing fields if present
    );
    boost::array<FlagDef,Algorithm::N_FLAGS> const & flagDefs = Algorithm::getFlagDefinitions();
    for (std::size_t i = 0; i < Algorithm::N_FLAGS; ++i) {
        _flags[i+1] = schema.addField(
            afw::table::Field<afw::table::Flag>(flagDefs[i].name, flagDefs[i].doc),
            true // replace existing fields if present
        );
    }
}

template <typename Algorithm>
void FlagsComponentMapper<Algorithm>::fail(
    afw::table::BaseRecord & record,
    MeasurementError const * error
) const {
    if (error) {
        assert(error->getFlagBit() < Algorithm::N_FLAGS);
        record.set(_flags[error->getFlagBit() + 1], true);
    }
    record.set(_flags[0], true);
}

template <typename Algorithm>
void FlagsComponentMapper<Algorithm>::apply(
    afw::table::BaseRecord & record,
    FlagsComponent<Algorithm> const & result
) const {
    for (std::size_t i = 0; i < Algorithm::N_FLAGS; ++i) {
        record.set(_flags[i+1], result._flags[i]);
    }
    record.set(_flags[0], false);
}

// Object that transfers values from FluxComponent to a record object.
class FluxComponentMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    FluxComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, FluxComponent const & result) const;

private:
    afw::table::Key<Flux> _flux;
    afw::table::Key<ErrElement> _fluxSigma;
};

// Object that transfers values from CentroidComponent to a record object.
class CentroidComponentMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    CentroidComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, CentroidComponent const & result) const;

private:
    afw::table::Key<CentroidElement> _x;
    afw::table::Key<CentroidElement> _y;
    afw::table::Key<ErrElement> _xSigma;
    afw::table::Key<ErrElement> _ySigma;
    afw::table::Key<ErrElement> _x_y_Cov;
};

// Object that transfers values from ShapeComponent to a record object.
class ShapeComponentMapper {
public:

    // Construct the mapper, creating fields by prepending the prefix to a set of standard names.
    ShapeComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, ShapeComponent const & result) const;

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
struct ResultMapper1 : public T1, public FlagsComponentMapper<Algorithm> {

    ResultMapper1(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1
    ) : T1(schema, prefix, uncertainty1),
        FlagsComponentMapper<Algorithm>(schema, prefix)
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        FlagsComponentMapper<Algorithm>::apply(record, result);
    }

};

template <typename Algorithm, typename T1, typename T2>
struct ResultMapper2 : public T1, public T2, public FlagsComponentMapper<Algorithm> {

    ResultMapper2(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1,
        ResultMapperUncertaintyEnum uncertainty2
    ) : T1(schema, prefix, uncertainty1),
        T2(schema, prefix, uncertainty2),
        FlagsComponentMapper<Algorithm>(schema, prefix)
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        T2::apply(record, result);
        FlagsComponentMapper<Algorithm>::apply(record, result);
    }

};

template <typename Algorithm, typename T1, typename T2, typename T3>
struct ResultMapper3 : public T1, public T2, public T3, public FlagsComponentMapper<Algorithm> {

    ResultMapper3(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1,
        ResultMapperUncertaintyEnum uncertainty2,
        ResultMapperUncertaintyEnum uncertainty3
    ) : T1(schema, prefix, uncertainty1),
        T2(schema, prefix, uncertainty2),
        T3(schema, prefix, uncertainty3),
        FlagsComponentMapper<Algorithm>(schema, prefix)
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        T2::apply(record, result);
        T3::apply(record, result);
        FlagsComponentMapper<Algorithm>::apply(record, result);
    }

};

template <typename Algorithm, typename T1, typename T2, typename T3, typename T4>
struct ResultMapper4 : public T1, public T2, public T3, public T4, public FlagsComponentMapper<Algorithm> {

    ResultMapper4(
        afw::table::Schema & schema,
        std::string const & prefix,
        ResultMapperUncertaintyEnum uncertainty1,
        ResultMapperUncertaintyEnum uncertainty2,
        ResultMapperUncertaintyEnum uncertainty3,
        ResultMapperUncertaintyEnum uncertainty4
    ) : T1(schema, prefix, uncertainty1),
        T2(schema, prefix, uncertainty2),
        T3(schema, prefix, uncertainty3),
        T4(schema, prefix, uncertainty4),
        FlagsComponentMapper<Algorithm>(schema, prefix)
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        T2::apply(record, result);
        T3::apply(record, result);
        T4::apply(record, result);
        FlagsComponentMapper<Algorithm>::apply(record, result);
    }

};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ResultMappers_h_INCLUDED
