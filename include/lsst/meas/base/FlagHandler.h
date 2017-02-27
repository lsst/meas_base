/*// -*- lsst-c++ -*-
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
 *  @brief Simple class used to define and document flags
 *         The name and doc constitute the identity of the FlagDefinition
 *         The number is used for indexing, but is assigned arbitrarily
*/
struct FlagDefinition {

    static const std::size_t number_undefined = SIZE_MAX;

    FlagDefinition() {
    }

    FlagDefinition(std::string _name, std::string _doc) {
        name = _name;
        doc = _doc;
        number = number_undefined;
    }

    FlagDefinition(std::string _name, std::string _doc, std::size_t  _number) {
        name = _name;
        doc = _doc;
        number = _number;
    }

    // equality of this type is based solely on the name key
    bool operator==(FlagDefinition const & other) const {
        return (other.name == name);
    }

    std::string name;
    std::string doc;
    std::size_t number;
};

/**
 *  @brief vector-type utility class to build a collection of FlagDefinitions
 */
class FlagDefinitionList {
public:
    /**
     *  @brief initialize a FlagDefinition list with no entries.
     */
    FlagDefinitionList() {
    };

#ifndef SWIG
    /**
     *  @brief initialize a FlagDefinition list from initializer_list.
     */
    FlagDefinitionList(std::initializer_list<FlagDefinition> const & list) {
        for (FlagDefinition const * iter = list.begin(); iter < list.end(); iter++) {
            add(iter->name, iter->doc);
        }
    }
#endif

    static FlagDefinitionList const & getEmptyList() {
        static FlagDefinitionList list;
        return list;
    }
    /**
     *  @brief get a reference to the FlagDefinition with specified index.
     */
    FlagDefinition getDefinition(size_t index) const {
        return _vector[index];
    }
    /**
     *  @brief get a reference to the FlagDefinition with specified array index
     */
    FlagDefinition operator[](size_t index) const {
        return getDefinition(index);
    }
    /**
     *  @brief get a reference to the FlagDefinition with specified name.
     */
    FlagDefinition getDefinition(std::string const & name) const {
        for (std::size_t i = 0; i < size(); i++) {
            if (_vector[i].name == name) return _vector[i];
        }
        throw FatalAlgorithmError("No Flag Definition for " + name);
    }
    /**
     *  @brief See if there is a FlagDefinition with specified name.
     */
    bool hasDefinition(std::string const & name) const {
        for (std::size_t i = 0; i < size(); i++) {
            if (_vector[i].name == name) return true;
        }
        return false;
    }
    /**
     *  @brief Add a Flag Defintion to act as a "General" failure flag
     *  This flag will be set if a Measurement error is thrown
     */
    FlagDefinition addFailureFlag(std::string const & doc="General Failure Flag");

    /**
     *  @brief Add a new FlagDefinition to this list. Return a copy with the
     *  FlagDefinition.number set corresponding to its index in the list.
     */
    FlagDefinition add(std::string const name, std::string const doc) {
        FlagDefinition flagDef = FlagDefinition(name, doc, _vector.size());
        _vector.push_back(flagDef);
        return _vector.back();
    }
    /**
     *  @brief return the current size (number of defined elements) of the collection
     */

    std::size_t size() const { return _vector.size(); }

private:
    mutable std::vector<FlagDefinition> _vector;
};

/**
 *  Utility class for handling flag fields that indicate the failure modes of an algorithm.
 *
 *  The typical pattern for using FlagHandler within an Algorithm is:
 *
 *   - Add a FlagHandler object as a data member.
 *
 *   - Create a FlagDefinitionList to specify the name and doc for each flag
 *     Add a "general" failure flag if one is needed (this indicates that some failure occurred).
 *     Add specific error flags for each type of error (these indicate a specific failure).
 *
 *   - Initialize the FlagHandler data member within the Algorithm's constructor,
 *     using the static addFields method to add the flags from the FlagDefinitionList
 *     to the schema.
 *
 *  See PsfFluxAlgorithm for a complete example.
 */
class FlagHandler {
public:
    /**
     *  Each error should have a corresponding static FlagDefinition object.
     *  In the Algorithm header file, this will be defined like this:
     *
     *      static FlagDefinition const & FAILURE;
     *      static FlagDefinition const & SOME_OTHER_FAILURE_MODE;
     *          ...
     *
     *  A static FlagDefinitionList is created in the Algorithm .cc file, like this:
     *
     *      FlagDefinitionList flagDefinitions;
     *      FlagDefinition const FAILURE = flagDefinitions.addFailureFlag();
     *      FlagDefinition const FAILURE_MODE = flagDefinitions.add("flag_mode", "Specific failure flag");
     *  @endcode
     */

    /**
     *  Default constructor for delayed initialization.
     *
     *  This constructor creates an invalid, unusable FlagHandler in the same way an const_iterator
     *  default constructor constructs an invalid const_iterator.  Its only purpose is to delay construction
     *  of the FlagHandler from an Algorithm constructor's initializer list to the constructor body,
     *  which can be necessary when the list of possible flags depends on the algorithm's configuration.
     *  To use this constructor to delay initialization, simply use it in the initializer list, and then
     *  assign the result of a call to addFields() to the FlagHandler data member later in the constructor.
     */
    FlagHandler() : failureFlagNumber(FlagDefinition::number_undefined) {}

    /**
     *  Define the universal name of the general failure flag
    */
    static std::string const & getFailureFlagName() {
        static std::string name = "flag";
        return name;
    }
    /**
     *  Add Flag fields to a schema, creating a FlagHandler object to manage them.
     *
     *  This is the way FlagHandlers will typically be constructed for new algorithms.
     *
     *  @param[out]  schema    Schema to which fields should be added.
     *  @param[in]   prefix    String name of the algorithm or algorithm component.  Field names will
     *                         be constructed by using schema.join() on this and the flag name from the
     *                         FlagDefinition array.
     *  @param[in]   flagDefs  Reference to a FlagDefinitionList
     *  @param[in]   exclDefs  optional FlagDefinitionList of flags to exclude
     *
     *  If the set of flags depends on the algorithm configuration, a flag may be excluded from the
     *  schema using the optional exclDefs parameter.  This can be specified using an initializer_list,
     *  as in: _flagHandler = FlagHandler::addFields(schema, prefix, flagDefs, {NO_PSF})
     */
    static FlagHandler addFields(
        afw::table::Schema & schema,
        std::string const & prefix,
        FlagDefinitionList const & flagDefs,
        FlagDefinitionList const & exclDefs=FlagDefinitionList::getEmptyList()
    );
    /**
     *  Construct a FlagHandler to manage fields already added to a schema.
     *
     *  This is primarily intended for use by forced measurement algorithms that need to parse the flags
     *  of the single-frame measurement algorithms providing their reference parameters.
     *
     *  @param[in]   A SubSchema object that holds the fields to extract and their namespace.
     *               Obtainable from the arguments to addFields() as "schema[prefix]".
     *  @param[in]   flagDefs  Reference to a FlagDefinitionList
     *  @param[in]   exclDefs  optional FlagDefinitionList of flags to exclude
     *
     *  As with addFields(), pointers must be valid only for the duration of this constructor call.
     */
    FlagHandler(
        afw::table::SubSchema const & s,
        FlagDefinitionList const & flagDefs,
        FlagDefinitionList const & exclDefs=FlagDefinitionList::getEmptyList()
    );
    /**
     *  Return the index of a flag with the given flag name
    */
    unsigned int getFlagNumber(std::string const & flagName) const {
        for (unsigned int i=0; i < _vector.size(); i++) {
            if (_vector[i].first == flagName && _vector[i].second.isValid()) {
                return i;
            }
        }
        throw FatalAlgorithmError("No FlagHandler entry for " + flagName);
    }
    /**
     *  Return the value of the flag name corresponding to the given flag index.
     */
    std::string getFlagName(std::size_t i) const {
        if (i < _vector.size() && _vector[i].second.isValid()) {
            return _vector[i].first;
        }
        throw FatalAlgorithmError("No legal FlagHandler entry number " + std::to_string(i));
    }
    /**
     *  Return the value of the flag field corresponding to the given flag index.
     */
    bool getValue(afw::table::BaseRecord const & record, std::size_t i) const {
        if (i < _vector.size() && _vector[i].second.isValid()) {
            return record.get(_vector[i].second);
        }
        throw FatalAlgorithmError("No legal FlagHandler entry number " + std::to_string(i));
    }
    /**
     *  Return the value of the flag field with the given flag name
     */
    bool getValue(afw::table::BaseRecord const & record, std::string flagName) const {
        for (std::size_t i = 0; i < _vector.size(); i++) {
            if (_vector[i].first == flagName && _vector[i].second.isValid()) {
                return record.get(_vector[i].second);
            }
        }
        throw FatalAlgorithmError("No FlagHandler entry for " + flagName);
    }
    /**
     *  Set the flag field corresponding to the given flag index.
     */
    void setValue(afw::table::BaseRecord & record, std::size_t i, bool value) const {
        if (i < _vector.size() && _vector[i].second.isValid()) {
            record.set(_vector[i].second, value);
            return;
        }
        throw FatalAlgorithmError("No legal FlagHandler entry number " + std::to_string(i));
    }
    /**
     *  Set the flag field corresponding to the given FlagDefinition Address.
     */
    /**
     *  Set the flag field corresponding to the given flag name.
     */
    void setValue(afw::table::BaseRecord & record, std::string flagName, bool value) const {
        for (std::size_t i = 0; i < _vector.size(); i++) {
            if (_vector[i].first == flagName && _vector[i].second.isValid()) {
                record.set(_vector[i].second, value);
                return;
            }
        }
        throw FatalAlgorithmError("No FlagHandler entry for " + flagName);
    }
    /**
     *  Get the index of the General Failure flag, if one is defined.  This flag is defined
     *  by most algorithms, and if defined, is set whenever an error is caught by the FlagHandler.
     *  If no General Failure flag is defined, this routine will return FlagDefinition::number_undefined 
     */
    std::size_t getFailureFlagNumber() const {
        return failureFlagNumber;
    }
    /**
     *  Handle an expected or unexpected Exception thrown by a measurement algorithm.
     *
     *  If the exception is expected, it should inherit from MeasurementError and can be passed here;
     *  this allows handleFailure to extract the failure mode enum value from the exception and set
     *  the corresponding flag.  The general failure flag will be set regardless of whether the "error"
     *  argument is NULL (which happens when an unexpected error occurs).
     */
    void handleFailure(afw::table::BaseRecord & record, MeasurementError const * error=NULL) const;

    std::size_t failureFlagNumber;
private:
    typedef std::vector< std::pair<std::string, afw::table::Key<afw::table::Flag> > > Vector;
    Vector _vector;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_FlagHandler_h_INCLUDED
