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

#ifndef LSST_MEAS_BASE_SdssShape_h_INCLUDED
#define LSST_MEAS_BASE_SdssShape_h_INCLUDED

#include <bitset>

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Algorithm.h"

namespace lsst { namespace meas { namespace base {

class SdssShapeResult;

/**
 *  @brief A C++ control class to handle SdssShapeAlgorithm's configuration
 *
 *  @copydetails PsfFluxControl
 */
class SdssShapeControl {
public:
    LSST_CONTROL_FIELD(background, double, "Additional value to add to background");
    LSST_CONTROL_FIELD(maxIter, int, "Maximum number of iterations");
    LSST_CONTROL_FIELD(maxShift, double, "Maximum centroid shift, limited to 2-10");
    LSST_CONTROL_FIELD(tol1, float, "Convergence tolerance for e1,e2");
    LSST_CONTROL_FIELD(tol2, float, "Convergence tolerance for FWHM");

    /// @copydoc PsfFluxControl::PsfFluxControl
    SdssShapeControl() : background(0.0), maxIter(100), maxShift(), tol1(1E-5), tol2(1E-4) {}

};

/**
 *  @brief A FunctorKey that maps SdssShapeResult to afw::table Records.
 *
 *  This is used internally by SdssShapeAlgorithm to transfer results from SdssShapeResult to SourceRecord,
 *  but it can also be used in the other direction by codes that need to extra an SdssShapeResult from
 *  a record.
 */
class SdssShapeResultKey : public afw::table::FunctorKey<SdssShapeResult> {
public:

    /**
     *  @brief Add the appropriate fields to a Schema, and return a SdssShapeResultKey that manages them
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_xx", "_yy", etc. will be appended to this
     *                         to form the full field names.
     */
    static SdssShapeResultKey addFields(
        afw::table::Schema & schema,
        std::string const & name
    );

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    SdssShapeResultKey() {}

    /**
     *  @brief Construct from a subschema, assuming _xx, _yy, etc. subfields
     *
     *  If a schema has "a_xx", "a_yy", etc. fields, this constructor allows you to construct
     *  a SdssShapeResultKey via:
     *  @code
     *  SdssShapeResultKey k(schema["a"]);
     *  @endcode
     */
    SdssShapeResultKey(afw::table::SubSchema const & s);

    /// Get a CentroidResult from the given record
    virtual SdssShapeResult get(afw::table::BaseRecord const & record) const;

    /// Set a CentroidResult in the given record
    virtual void set(afw::table::BaseRecord & record, SdssShapeResult const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Keys
    bool operator==(SdssShapeResultKey const & other) const;
    bool operator!=(SdssShapeResultKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if the key is valid.
    bool isValid() const;

    FlagHandler const & getFlagHandler() const { return _flagHandler; }

private:
    ShapeResultKey _shapeResult;
    CentroidResultKey _centroidResult;
    FluxResultKey _fluxResult;
    afw::table::Key<ShapeElement> _xy4;
    afw::table::Key<ErrElement> _xy4Sigma;
    afw::table::Key<ErrElement> _flux_xx_Cov;
    afw::table::Key<ErrElement> _flux_yy_Cov;
    afw::table::Key<ErrElement> _flux_xy_Cov;
    FlagHandler _flagHandler;
};

/**
 *  @brief Measure the image moments of source using adaptive Gaussian weights.
 *
 *  This algorithm measures the weighted second moments of an image using a Gaussian weight function, which
 *  is iteratively updated to match the current weights.  If this iteration does not converge, it can fall
 *  back to using unweighted moments, which can be significantly noisier.
 */
class SdssShapeAlgorithm : public SimpleAlgorithm {
public:

    typedef SdssShapeControl Control;
    typedef SdssShapeResult Result;
    typedef SdssShapeResultKey ResultKey;

    enum {
        FAILURE=FlagHandler::FAILURE,
        UNWEIGHTED_BAD,
        UNWEIGHTED,
        SHIFT,
        MAXITER,
        N_FLAGS
    };

    SdssShapeAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    template <typename T>
    static Result apply(
        afw::image::MaskedImage<T> const & image,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    template <typename T>
    static Result apply(
        afw::image::Image<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

private:

    Control _ctrl;
    ResultKey _resultKey;
    SafeCentroidExtractor _centroidExtractor;
};

/**
 *  @brief Result object SdssShapeAlgorithm
 *
 *  Because we have use cases for running SdssShape outside of the measurement framework (in particular,
 *  we need to run it on PSF model images), we provide an interface that doesn't need to use SourceRecord
 *  for its inputs and outputs.  Instead, it returns an instance of this class.
 *
 *  Note: for what I guess are historical reasons, SdssShape computes covariance terms between the flux
 *  and the shape, but not between the flux and centroid or centroid and shape.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class SdssShapeResult : public ShapeResult, public CentroidResult, public FluxResult {
public:
    ShapeElement xy4;       ///< A fourth moment used in lensing (RHL needs to clarify; not in the old docs)
    ErrElement xy4Sigma;    ///< 1-Sigma uncertainty on xy4
    ErrElement flux_xx_Cov; ///< flux, xx term in the uncertainty covariance matrix
    ErrElement flux_yy_Cov; ///< flux, yy term in the uncertainty covariance matrix
    ErrElement flux_xy_Cov; ///< flux, xy term in the uncertainty covariance matrix

#ifndef SWIG
    std::bitset<SdssShapeAlgorithm::N_FLAGS> flags; ///< Status flags (see SdssShapeAlgorithm).
#endif

    /// Flag getter for Swig, which doesn't understand std::bitset
    bool getFlag(int index) const { return flags[index]; }

    SdssShapeResult(); ///< Constructor; initializes everything to NaN

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssShape_h_INCLUDED
