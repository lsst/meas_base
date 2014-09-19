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

#ifndef LSST_MEAS_BASE_exceptions_h_INCLUDED
#define LSST_MEAS_BASE_exceptions_h_INCLUDED

#include "lsst/pex/exceptions.h"

namespace lsst { namespace meas { namespace base {



/**
 *  @brief Exception to be thrown when a measurement algorithm experiences a known failure mode.
 *
 *  In addition to the usual message, MeasurementError must be constructed with the bit of the
 *  algorithm-specific flag that indicates the known failure mode.  This allows the measurement
 *  framework to set that flag upon failure.  Typically, this flag bit is also used to look up
 *  the message from the algorithm classes FlagDefinition list; the common pattern is:
 *  @code
 *      if (badThingHappened) {
 *          // BAD_THING is an enum value from the Algorithm's FlagBits enum; getFlagDefinitions()
 *          // is a static method algorithms are expected to define.
 *          throw LSST_EXCEPT(MeasurementError, getFlagDefinitions()[BAD_THING), BAD_THING);
 *       }
 *  @endcode
 */
class MeasurementError : public pex::exceptions::RuntimeError {
public:

    /// Constructor; should only be invoked from Python macro
    MeasurementError(std::string const & message, std::size_t flagBit) :
        pex::exceptions::RuntimeError(message),
        _flagBit(flagBit)
    {}

    /// Constructor; should only be invoked by the LSST_EXCEPT macro (see class docs)
    MeasurementError(LSST_EARGS_TYPED, std::size_t flagBit) :
        pex::exceptions::RuntimeError(LSST_EARGS_UNTYPED),
        _flagBit(flagBit)
    {}

    /// Return the flag bit associated with the error.
    std::size_t getFlagBit() const { return _flagBit; }

    virtual char const* getType(void) const throw() { return "lsst::meas::base::MeasurementError *"; };

    virtual lsst::pex::exceptions::Exception* clone(void) const {
        return new MeasurementError(*this);
    };

private:
    std::size_t _flagBit;
};

/**
 *  @brief Exception to be thrown when a measurement algorithm experiences a fatal error.
 *
 *  This error type causes the meas_base framework to throw completely out of the measurement loop
 *  which is run for each exposure, sourceCatalog pair.
 */
class FatalAlgorithmError : public pex::exceptions::RuntimeError {
public:

    /// Constructor; should only be invoked from Python macro
    FatalAlgorithmError(std::string const & message) :
        pex::exceptions::RuntimeError(message)
    {}

    /// Constructor; should only be invoked by the LSST_EXCEPT macro (see class docs)
    FatalAlgorithmError(LSST_EARGS_TYPED) :
        pex::exceptions::RuntimeError(LSST_EARGS_UNTYPED)
    {}

    virtual char const* getType(void) const throw() { return "lsst::meas::base::FatalAlgorithmError *"; };

    virtual lsst::pex::exceptions::Exception* clone(void) const {
        return new FatalAlgorithmError(*this);
    };
};

/**
 *  @brief Exception to be thrown when a measurement algorithm encounters a NaN or infinite pixel.
 *
 *  When caught by the plugin framework, this exception will generate a log message.
 */
LSST_EXCEPTION_TYPE(PixelValueError, lsst::pex::exceptions::DomainError,
                    lsst::meas::base::PixelValueError);

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_exceptions_h_INCLUDED
