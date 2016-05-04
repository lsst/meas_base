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

/**
 *  Utility class for handling flag fields that indicate the failure modes of an algorithm.
 *
 *  The typical pattern for using FlagHandler within an Algorithm is:
 *   - Add a FlagHandler object as a data member.
 *   - Create an old-style enum that defines all of the failure modes.  This must start with a
 *     "general" failure flag that is set on any fatal condition, set to the numeric value
 *     FlagHandler::FAILURE.
 *   - Create a static array of FlagDefinition to hold the field names for the error flags, with each
 *     entry corresponding to an entry in the enum.  Names here should not include the algorihtm name,
 *     and should generally have the form "flag_*", with the first entry (corresponding to the general
 *     failure flag) simply "flag".
 *   - Initialize the FlagHandler data member within the Algorithm's constructor, using the addFields
 *     static method to add fields to the schema at the same time the FlagHandler is constructed to
 *     manage them.
 *
 *  See PsfFluxAlgorithm for a complete example.
 */
class FlagHandler {
public:

    /**
     *  Required enum values for all Algorithm failure mode enumerations.
     *
     *  This should be considered logically to be a "base class" for all Algorithm failure enums;
     *  enums can't actually have a base class in C/C++, but we can fake this with the following pattern:
     *  @code
     *  class MyAlgorithm : public Algorithm {
     *  public:
     *      enum {
     *          FAILURE=FlagHandler::FAILURE,
     *          SOME_OTHER_FAILURE_MODE,
     *          ...
     *      };
     *      ...
     *  };
     *  @endcode
     */
    enum { FAILURE=0 };

    /**
     *  Default constructor for delayed initialization.
     *
     *  This constructor creates an invalid, unusable FlagHandler in the same way an iterator
     *  default constructor constructs an invalid iterator.  Its only purpose is to delay construction
     *  of the FlagHandler from an Algorithm constructor's initializer list to the constructor body,
     *  which can be necessary when the list of possible flags depends on the algorithm's configuration.
     *  To use this constructor to delay initialization, simply use it in the initializer list, and then
     *  assign the result of a call to addFields() to the FlagHandler data member later in the constructor.
     */
    FlagHandler() {}

    /**
     *  Add Flag fields to a schema, creating a FlagHandler object to manage them.
     *
     *  This is the way FlagHandlers will typically be constructed for new algorithms.
     *
     *  @param[out]  schema    Schema to which fields should be added.
     *  @param[in]   prefix    String name of the algorithm or algorithm component.  Field names will
     *                         be constructed by using schema.join() on this and the flag name from the
     *                         FlagDefinition array.
     *  @param[in]   begin     Iterator to the beginning of an array of FlagDefinition.
     *  @param[in]   end       Iterator to one past the end of an array of FlagDefinition.
     *
     *  We use pointers rather than an iterator type for the FlagDefinition array to allow the user
     *  maximum flexibility in the array type - C arrays, std::array, std::array, and std::vector
     *  (as well as any other container with contiguous memory) may be used.  The pointers must remain
     *  valid only for the duration of the call to this function.
     */
    static FlagHandler addFields(
        afw::table::Schema & schema,
        std::string const & prefix,
        FlagDefinition const * begin,
        FlagDefinition const * end
    );

    /**
     *  Construct a FlagHandler to manage fields already added to a schema.
     *
     *  This is primarily intended for use by forced measurement algorithms that need to parse the flags
     *  of the single-frame measurement algorithms providing their reference parameters.
     *
     *  @param[in] s      A SubSchema object that holds the fields to extract and their namespace.
     *                    Obtainable from the arguments to addFields() as "schema[prefix]".
     *  @param[in] begin  Iterator to the beginning of an array of FlagDefinition.
     *  @param[in] end    Iterator to one past the end of an array of FlagDefinition.
     *
     *  As with addFields(), pointers must be valid only for the duration of this constructor call.
     */
    FlagHandler(
        afw::table::SubSchema const & s,
        FlagDefinition const * begin,
        FlagDefinition const * end
    );

    /**
     *  Return the FlagDefinition object that corresponds to a particular enum value.
     */
    FlagDefinition getDefinition(std::size_t i) const {
        assert(_vector.size() > i);  // Flag 'i' needs to be available
        return _vector[i].first;
    }

    /**
     *  Return the value of the flag field corresponding to the given enum value.
     */
    bool getValue(afw::table::BaseRecord const & record, std::size_t i) const {
        assert(_vector.size() > i);  // Flag 'i' needs to be available
        return record.get(_vector[i].second);
    }

    /**
     *  Set the flag field corresponding to the given enum value.
     */
    void setValue(afw::table::BaseRecord & record, std::size_t i, bool value) const {
        assert(_vector.size() > i);  // Flag 'i' needs to be available
        record.set(_vector[i].second, value);
    }

    /**
     *  Handle an expected or unexpected Exception thrown by a measurement algorithm.
     *
     *  If the exception is expected, it should inherit from MeasurementError and can be passed here;
     *  this allows handleFailure to extract the failure mode enum value from the exception and set
     *  the corresponding flag.  The general failure flag will be set regardless of whether the "error"
     *  argument is null.
     */
    void handleFailure(afw::table::BaseRecord & record, MeasurementError const * error=NULL) const;

private:

    typedef std::vector< std::pair<FlagDefinition, afw::table::Key<afw::table::Flag> > > Vector;

    Vector _vector;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_FlagHandler_h_INCLUDED
