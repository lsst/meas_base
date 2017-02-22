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
 *  @brief Simple POD struct used to define and document flags
 *         The name and doc constitute the identity of the FlagDefinition
 *         The number is used for indexing, but is assigned arbitrarily 
*/
struct FlagDefinition {

    static const std::size_t number_undefined = 1000000;
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

    bool operator==(FlagDefinition const & other) const {
        return (other.name == name && other.doc == doc);
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
     *  @brief initialize a FlagDefinition collection with no entries.
     */
    FlagDefinitionList() {
    };

#ifndef SWIG

    /**
     *  @brief initialize a FlagDefinition collection from initializer_list.
     */
    FlagDefinitionList(std::initializer_list<FlagDefinition> const & list) {
        for (FlagDefinition const * iter = list.begin(); iter < list.end(); iter++) {
            add(iter->name, iter->doc);
        }
    }
#endif

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
        throw lsst::pex::exceptions::RuntimeError("No Flag Definition for " + name);
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
     *  @brief Add a new FlagDefinition to this collection.
     *  Return a reference to the newly create definition with the
     *  FlagDefinition.number corresponding to its index in the array..
     */
    FlagDefinition add(std::string const name, std::string const doc) {
        FlagDefinition * flagDef = new FlagDefinition(name, doc, _vector.size());
        _vector.push_back(* flagDef);
        return _vector.back(); //* flagDef; 
    }
    FlagDefinition add(FlagDefinition const & flagDef) { return add(flagDef.name, flagDef.doc); }
    /**
     *  @brief return the current size (number of defined elements) of the collection
     */

    static FlagDefinition getFailureFlag() {
        static FlagDefinition flagDef = FlagDefinition("flag", "General Failure");
        return flagDef;
    }

    FlagDefinition addFailureFlag() {
        return add(getFailureFlag());
    }

    std::size_t size() const { return _vector.size(); }

private:
    mutable std::vector<FlagDefinition> _vector;
};

/**
 *  Utility class for handling flag fields that indicate the failure modes of an algorithm.
 *
 *  The typical pattern for using FlagHandler within an Algorithm is:
 *   - Add a FlagHandler object as a data member.
 *   - Create a FlagDefinitionList to hold the names and docs for the error flags, with each
 *     entry corresponding to   Names here should not include the algorihtm name,
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
     *  Each error should have a corresponding static FlagDefinition object.
     *  In the Algorithm header file, this will be defined like this:
     *
     *      static FlagDefinition const & FAILURE;
     *      static FlagDefinition const & SOME_OTHER_FAILURE_MODE;
     *          ...
     *
     *  A FlagDefinitionList is created in the .cc file and is used to create the FlagDefinition
     *  objects and insure that they are properly sequence:
     *
     *      FlagDefinitionList flagDefinitions;
     *      FlagDefinition const FAILURE = flagDefinitions.add(FlagHandler::FAILURE); 
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
     *  @param[in]   flagDefs  Reference to a FlagDefinitionList
     *  @param[in]   skipDefs  Reference to a FlagDefinitionList of entries in flagDefs skip
     *
     */
    static FlagHandler addFields(
        afw::table::Schema & schema,
        std::string const & prefix,
        FlagDefinitionList const & flagDefs,
        FlagDefinitionList const & skipDefs=FlagDefinitionList()
    );
    /**
     *  Construct a FlagHandler to manage fields already added to a schema.
     *
     *  This is primarily intended for use by forced measurement algorithms that need to parse the flags
     *  of the single-frame measurement algorithms providing their reference parameters.
     *
     *  @param[in] s      A SubSchema object that holds the fields to extract and their namespace.
     *                    Obtainable from the arguments to addFields() as "schema[prefix]".
     *  @param[in]   flagDefs  Reference to a FlagDefinitionList
     *  @param[in]   skipDefs  Reference to a FlagDefinitionList of entries in flagDefs skip
     *
     *  As with addFields(), pointers must be valid only for the duration of this constructor call.
     */
    FlagHandler(
        afw::table::SubSchema const & s,
        FlagDefinitionList const & flagDefs,
        FlagDefinitionList const & skipDefs=FlagDefinitionList()
    );
    /**
     *  Return the FlagDefinition object that corresponds to given FlagHandler name
     */
    FlagDefinition getDefinition(std::string const & flagName) const {
        for (std::size_t i = 0; i < _vector.size(); i++) {
            if (_vector[i].first.name == flagName && _vector[i].second.isValid()) {
                return _vector[i].first;
            }
        }
        throw pex::exceptions::RuntimeError("No FlagHandler entry for " + flagName);
    }
    /**
     *  Return the FlagDefinition object that corresponds to given FlagHandler index
     */
    FlagDefinition getDefinition(std::size_t i) const {
        if (i < _vector.size() && _vector[i].second.isValid()) {
            return _vector[i].first;
        }
        throw pex::exceptions::RuntimeError("No legal FlagHandler entry number " + std::to_string(i));
    }
    /**
     *  Return the index of a flag with the given flag name
    */
    unsigned int getFlagNumber(std::string const & flagName) const {
        for (unsigned int i=0; i < _vector.size(); i++) {
            if (_vector[i].first.name == flagName && _vector[i].second.isValid()) {
                return i;
            }
        }
        throw pex::exceptions::RuntimeError("No FlagHandler entry for " + flagName);
    }
    /**
     *  Return the value of the flag field corresponding to the given flag index.
     */
    bool getValue(afw::table::BaseRecord const & record, std::size_t i) const {
        if (i < _vector.size() && _vector[i].second.isValid()) {
            return record.get(_vector[i].second);
        }
        throw pex::exceptions::RuntimeError("No legal FlagHandler entry number " + std::to_string(i));
    }
    /**
     *  Return the value of the flag field with the given flag name
     */
    bool getValue(afw::table::BaseRecord const & record, std::string flagName) const {
        for (std::size_t i = 0; i < _vector.size(); i++) {
            if (_vector[i].first.name == flagName && _vector[i].second.isValid()) {
                return record.get(_vector[i].second);
            }
        }
        throw pex::exceptions::RuntimeError("No FlagHandler entry for " + flagName);
    }
    /**
     *  Set the flag field corresponding to the given flag index.
     */
    void setValue(afw::table::BaseRecord & record, std::size_t i, bool value) const {
        if (i < _vector.size() && _vector[i].second.isValid()) {
            record.set(_vector[i].second, value);
            return;
        }
        throw pex::exceptions::RuntimeError("No legal FlagHandler entry number " + std::to_string(i));
    }
    /**
     *  Set the flag field corresponding to the given FlagDefinition Address.
     */
    /**
     *  Set the flag field corresponding to the given flag name.
     */
    void setValue(afw::table::BaseRecord & record, std::string flagName, bool value) const {
        for (std::size_t i = 0; i < _vector.size(); i++) {
            if (_vector[i].first.name == flagName && _vector[i].second.isValid()) {
                record.set(_vector[i].second, value);
                return;
            }
        }
        throw pex::exceptions::RuntimeError("No FlagHandler entry for " + flagName);
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
