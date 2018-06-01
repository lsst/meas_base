// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST.
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
#include "lsst/meas/base/Transform.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/table/FunctorKey.h"
#include "lsst/afw/table/Schema.h"

namespace lsst {
namespace meas {
namespace base {

/**
 *  @brief A reusable result struct for flux measurements.
 */
struct FluxResult {
    Flux flux;                 ///< Measured flux in DN.
    FluxErrElement fluxSigma;  ///< 1-Sigma error (sqrt of variance) on flux in DN.

    /// Default constructor; initializes everything to NaN.
    FluxResult();

    /// Constructor from flux and its uncertainty
    explicit FluxResult(Flux flux_, FluxErrElement fluxSigma_) : flux(flux_), fluxSigma(fluxSigma_) {}
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
     *  The unit for both fields will be "count".
     */
    static FluxResultKey addFields(afw::table::Schema& schema, std::string const& name,
                                   std::string const& doc);

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    FluxResultKey() : _flux(), _fluxSigma() {}

    /// Construct from a pair of Keys
    FluxResultKey(afw::table::Key<meas::base::Flux> const& flux,  // namespace qualification to unconfuse swig
                  afw::table::Key<FluxErrElement> const& fluxSigma)
            : _flux(flux), _fluxSigma(fluxSigma) {}

    /**
     *  @brief Construct from a subschema, assuming flux and fluxSigma subfields
     *
     *  If a schema has "a_flux" and "a_fluxSigma" fields, this constructor allows you to construct
     *  a FluxResultKey via:
     *  @code
     *  FluxResultKey k(schema["a"]);
     *  @endcode
     */
    FluxResultKey(afw::table::SubSchema const& s) : _flux(s["flux"]), _fluxSigma(s["fluxSigma"]) {}

    /// Get a FluxResult from the given record
    virtual FluxResult get(afw::table::BaseRecord const& record) const;

    /// Set a FluxResult in the given record
    virtual void set(afw::table::BaseRecord& record, FluxResult const& other) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying flux and fluxSigma Keys
    bool operator==(FluxResultKey const& other) const {
        return _flux == other._flux && _fluxSigma == other._fluxSigma;
    }
    bool operator!=(FluxResultKey const& other) const { return !(*this == other); }
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

/**
 *  @brief A reusable result struct for magnitudes
 */
struct MagResult {
    Mag mag;
    MagErrElement magErr;
};

/**
 *  @brief A FunctorKey for MagResult
 *
 *  This class makes it easy to copy magnitudes and their uncertainties to and from records, and provides
 *  a method to add the appropriate fields to a Schema.
 */
class MagResultKey : public afw::table::FunctorKey<MagResult> {
public:
    /**
     *  Add a pair of _mag, _magErr fields to a Schema, and return a MagResultKey that points to them.
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_mag", "_magErr" will be appended to this
     *                         to form the full field names.
     */
    static MagResultKey addFields(afw::table::Schema& schema, std::string const& name);

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    MagResultKey() : _magKey(), _magErrKey() {}

    /**
     *  @brief Construct from a subschema, assuming mag and magErr subfields
     *
     *  If a schema has "a_mag" and "a_magErr" fields, this enables construction of a MagResultKey via:
     *  @code
     *  MagResultKey k(schema["a"]);
     *  @endcode
     */
    MagResultKey(afw::table::SubSchema const& s) : _magKey(s["mag"]), _magErrKey(s["magErr"]) {}

    /// Get a MagResult from the given record.
    virtual MagResult get(afw::table::BaseRecord const& record) const;

    /// Set a MagResult in the given record.
    virtual void set(afw::table::BaseRecord& record, MagResult const& magResult) const;

    /// Set a MagResult in the record given the result of `afw::image::Calib::getMagnitude(double, double)`.
    virtual void set(afw::table::BaseRecord& record, std::pair<double, double> const& magPair) const;

private:
    afw::table::Key<Mag> _magKey;
    afw::table::Key<MagErrElement> _magErrKey;
};

/**
 *  Base for flux measurement transformations
 *
 *  Provides a basic transform from flux plus associated uncertainty to
 *  magnitude with uncertainty. The basic "flag" attribute for the measurement
 *  algorithm is propagated to the output
 *
 *  Subclasses should define a constructor which take a Control argument
 *  corresponding to the measurement algorithm being transformed and ensure
 *  that any other necessary flags are propagated.
 */
class FluxTransform : public BaseTransform {
public:
    FluxTransform(std::string const& name, afw::table::SchemaMapper& mapper);

    /*
     * @brief Perform transformation from inputCatalog to outputCatalog.
     *
     * @param[in]     inputCatalog   Source of data to be transformed
     * @param[in,out] outputCatalog  Container for transformed results
     * @param[in]     wcs            World coordinate system under which transformation will take place
     * @param[in]     calib          Photometric calibration under which transformation will take place
     * @throws        LengthError    Catalog sizes do not match
     */
    virtual void operator()(afw::table::SourceCatalog const& inputCatalog,
                            afw::table::BaseCatalog& outputCatalog, afw::geom::SkyWcs const& wcs,
                            afw::image::Calib const& calib) const;

private:
    MagResultKey _magKey;
};

/**
 *  Temporarily replace negative fluxes with NaNs
 *
 *  Instantiating a NoThrowOnNegativeFluxContext object will cause afw::image::Calib
 *  objects to return NaN, rather than throwing, when asked to convert a negative flux to
 *  a magnitude for the lifetime of the NoThrowOnNegativeFluxContext.
 */
class NoThrowOnNegativeFluxContext {
public:
    NoThrowOnNegativeFluxContext();
    ~NoThrowOnNegativeFluxContext();

private:
    bool _throwOnNegative;
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_FluxUtilities_h_INCLUDED
