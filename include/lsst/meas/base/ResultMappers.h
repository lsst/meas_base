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
#include "lsst/afw/table/arrays.h"
#include "lsst/meas/base/Results.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @addtogroup measBaseResults
 *  @{
 */

/**
 *  @brief An object that transfers values from FlagsComponent to afw::table::BaseRecord
 *
 *  This is implicitly included in all of @ref measBaseResultMapperTemplates, and will generally only be
 *  used (directly) by these.  It provides the fail() method to those templates.
 */
template <typename Algorithm>
struct FlagsComponentMapper {

    /**
     *  @brief Construct the mapper, adding fields to the given schema and saving their keys
     *
     *  Field names will have the form "<prefix>_flag_<FlagDef.doc>", where "<prefix>" is the prefix
     *  argument to this constructor, and "<FlagDef.doc>" is the string documentation for this flag
     *  obtained by a call to Algorithm::getFlagDefinitions().
     */
    FlagsComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix
    );

    /**
     *  @brief Record a failure for the algorithm in the given record.
     *
     *  Sets the general failure flag field, and if error is not null, the specific bit held by the given
     *  exception.
     *
     *  This method is called by the measurement framework when an algorithm throws an exception, but
     *  the error argument is only non-null when that exception is a MeasurementError.  Algorithms that
     *  wish to set flags that do not indicate complete failure should set them themselves directly in the
     *  result object, as throwing MeasurementError (and hence invoking fail()) does not allow other values
     *  to be returned.
     */
    void fail(afw::table::BaseRecord & record, MeasurementError const * error=NULL) const;

    /// Transfer values from the result struct to the record, and clear the general failure flag field.
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
            afw::table::Field<afw::table::Flag>((prefix + "_flag_") + flagDefs[i].name, flagDefs[i].doc),
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


/**
 *  @brief An object that transfers values from FluxComponent to afw::table::BaseRecord
 *
 *  This should be included in one of @ref measBaseResultMapperTemplates to correspond with using
 *  FluxComponent in the same position in one of @ref measBaseResultTemplates, and will otherwise
 *  not be used directly by users.
 */
class FluxComponentMapper {
public:

    /**
     *  @brief Construct the mapper, adding fields to the given schema and saving their keys
     *
     *  The given prefix will form the first part of all fields, and the uncertainty argument
     *  sets which uncertainty fields will be added to the schema and transferred during apply().
     */
    FluxComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        UncertaintyEnum uncertainty
    );

    /// Transfer values from the result struct to the record
    void apply(afw::table::BaseRecord & record, FluxComponent const & result) const;

private:
    afw::table::Key<Flux> _flux;
    afw::table::Key<FluxErrElement> _fluxSigma;
};


/**
 *  @brief An object that transfers values from CentroidComponent to afw::table::BaseRecord
 *
 *  This should be included in one of @ref measBaseResultMapperTemplates to correspond with using
 *  CentroidComponent in the same position in one of @ref measBaseResultTemplates, and will otherwise
 *  not be used directly by users.
 */
class CentroidComponentMapper {
public:

    /**
     *  @brief Construct the mapper, adding fields to the given schema and saving their keys
     *
     *  The given prefix will form the first part of all fields, and the uncertainty argument
     *  sets which uncertainty fields will be added to the schema and transferred during apply().
     */
    CentroidComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        UncertaintyEnum uncertainty
    );

    // Transfer values from the result struct to the record
    void apply(afw::table::BaseRecord & record, CentroidComponent const & result) const;

private:
    afw::table::Key<CentroidElement> _x;
    afw::table::Key<CentroidElement> _y;
    afw::table::Key<ErrElement> _xSigma;
    afw::table::Key<ErrElement> _ySigma;
    afw::table::Key<ErrElement> _x_y_Cov;
};


/**
 *  @brief An object that transfers values from CentroidComponent to afw::table::BaseRecord
 *
 *  This should be included in one of @ref measBaseResultMapperTemplates to correspond with using
 *  CentroidComponent in the same position in one of @ref measBaseResultTemplates, and will otherwise
 *  not be used directly by users.
 */
class ShapeComponentMapper {
public:

    /**
     *  @brief Construct the mapper, adding fields to the given schema and saving their keys
     *
     *  The given prefix will form the first part of all fields, and the uncertainty argument
     *  sets which uncertainty fields will be added to the schema and transferred during apply().
     */
    ShapeComponentMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        UncertaintyEnum uncertainty
    );

    /// Transfer values from the result struct to the record
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

/**
 *  @defgroup measBaseResultMapperTemplates the ResultMapperN templates
 *
 *  These templates aggregate ComponentMapper objects in the same way that the @ref measBaseResultTemplates
 *  aggregate the Components themselves.  Each Algorithm class should define a @c ResultMapper typedef
 *  to one of these templates, with argments that correspond to the arguments in the @c Result typedef.
 *
 *  When an Algorithm defines a custom Component, it must also define a custom ComponentMapper.  This class
 *  must meet the following requirements (the standard ComponentMapper classes can all serve as examples):
 *    - It must be copyable (simply because ResultMappers are returned by value).
 *    - It must have a constructor with the same signature as FluxComponentMapper::FluxComponentMapper,
 *      which should add fields to the Schema object and save the returned keys for later use.
 *    - It must have an apply() method that takes a non-const reference to afw::table::BaseRecord and
 *      a const reference to the corresponding Component object, which transfers values from the custom
 *      Component to the record.
 *  @{
 */

template <typename Algorithm>
struct ResultMapper0 : public FlagsComponentMapper<Algorithm> {

    template <typename A1>
    ResultMapper0(
        afw::table::Schema & schema,
        std::string const & prefix,
        A1 a1
    ) : FlagsComponentMapper<Algorithm>(schema, prefix)
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        FlagsComponentMapper<Algorithm>::apply(record, result);
    }

};

template <typename Algorithm, typename T1>
struct ResultMapper1 : public T1, public FlagsComponentMapper<Algorithm> {

    template <typename A1>
    ResultMapper1(
        afw::table::Schema & schema,
        std::string const & prefix,
        A1 a1
    ) : T1(schema, prefix, a1),
        FlagsComponentMapper<Algorithm>(schema, prefix)
    {}

    void apply(afw::table::BaseRecord & record, typename Algorithm::Result const & result) {
        T1::apply(record, result);
        FlagsComponentMapper<Algorithm>::apply(record, result);
    }

};

template <typename Algorithm, typename T1, typename T2>
struct ResultMapper2 : public T1, public T2, public FlagsComponentMapper<Algorithm> {

    template <typename A1, typename A2>
    ResultMapper2(
        afw::table::Schema & schema,
        std::string const & prefix,
        A1 a1, A2 a2
    ) : T1(schema, prefix, a1),
        T2(schema, prefix, a2),
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

    template <typename A1, typename A2, typename A3>
    ResultMapper3(
        afw::table::Schema & schema,
        std::string const & prefix,
        A1 a1, A2 a2, A3 a3
    ) : T1(schema, prefix, a1),
        T2(schema, prefix, a2),
        T3(schema, prefix, a3),
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

    template <typename A1, typename A2, typename A3, typename A4>
    ResultMapper4(
        afw::table::Schema & schema,
        std::string const & prefix,
        A1 a1, A2 a2, A3 a3, A4 a4
    ) : T1(schema, prefix, a1),
        T2(schema, prefix, a2),
        T3(schema, prefix, a3),
        T4(schema, prefix, a4),
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

/** @} */ // end of measBaseResultMapperTemplates group

/** @} */ // end of measBaseResults group

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ResultMappers_h_INCLUDED
