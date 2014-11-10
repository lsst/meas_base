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

#ifndef LSST_MEAS_BASE_FluxUtilities_h_INCLUDED
#define LSST_MEAS_BASE_FluxUtilities_h_INCLUDED

#include "lsst/meas/base/constants.h"
#include "lsst/afw/table/FunctorKey.h"
#include "lsst/afw/table/Schema.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A reusable result struct for flux measurements.
 */
struct FluxResult {
    Flux flux; ///< Measured flux in DN.
    FluxErrElement fluxSigma; ///< 1-Sigma error (sqrt of variance) on flux in DN.

    /// Default constructor; initializes everything to NaN.
    FluxResult();

    /// Constructor from flux and its uncertainty
    explicit FluxResult(Flux flux_, FluxErrElement fluxSigma_) :
        flux(flux_), fluxSigma(fluxSigma_)
    {}
};

/**
 *  @brief A FunctorKey for FluxResult
 *
 *  This class makes it easy to copy fluxes and their uncertainties to and from records, and provides
 *  a method to add the appropriate fields to a Schema.
 */
class FluxResultKey : public afw::table::FunctorKey<FluxResult> {
public:

    /**
     *  Add a pair of _flux, _fluxSigma fields to a Schema, and return a FluxResultKey that points to them.
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_flux", "_fluxSigma" will be appended to this
     *                         to form the full field names.
     *  @param[in]     doc     String used as the documentation for the fields.
     *
     *  The unit for both fields will be "dn".
     */
    static FluxResultKey addFields(
        afw::table::Schema & schema,
        std::string const & name,
        std::string const & doc,
        UncertaintyEnum uncertainty=SIGMA_ONLY

    );

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    FluxResultKey() : _flux(), _fluxSigma() {}

    /// Construct from a pair of Keys
    FluxResultKey(
        afw::table::Key<meas::base::Flux> const & flux,  // namespace qualification to unconfuse swig
        afw::table::Key<FluxErrElement> const & fluxSigma
    ) :
        _flux(flux), _fluxSigma(fluxSigma)
    {}

    /**
     *  @brief Construct from a subschema, assuming flux and fluxSigma subfields
     *
     *  If a schema has "a_flux" and "a_fluxSigma" fields, this constructor allows you to construct
     *  a FluxResultKey via:
     *  @code
     *  FluxResultKey k(schema["a"]);
     *  @endcode
     */
    FluxResultKey(afw::table::SubSchema const & s) : _flux(s["flux"]), _fluxSigma(s["fluxSigma"]) {}

    /// Get a FluxResult from the given record
    virtual FluxResult get(afw::table::BaseRecord const & record) const;

    /// Set a FluxResult in the given record
    virtual void set(afw::table::BaseRecord & record, FluxResult const & other) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying flux and fluxSigma Keys
    bool operator==(FluxResultKey const & other) const {
        return _flux == other._flux && _fluxSigma == other._fluxSigma;
    }
    bool operator!=(FluxResultKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if both the flux and fluxSigma Keys are valid.
    bool isValid() const { return _flux.isValid() && _fluxSigma.isValid(); }

    /// Return the underlying flux Key
    afw::table::Key<meas::base::Flux> getFlux() const { return _flux; }

    /// Return the underlying fluxSigma Key
    afw::table::Key<FluxErrElement> getFluxSigma() const { return _fluxSigma; }

private:
    afw::table::Key<Flux> _flux;
    afw::table::Key<FluxErrElement> _fluxSigma;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_FluxUtilities_h_INCLUDED
