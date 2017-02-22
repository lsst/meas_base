// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_BASE_SdssShape_h_INCLUDED
#define LSST_MEAS_BASE_SdssShape_h_INCLUDED

#include <bitset>

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/table/aggregates.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Algorithm.h"

namespace lsst {
namespace meas {
namespace base {

class SdssShapeResult;

/**
 *  @brief A C++ control class to handle SdssShapeAlgorithm's configuration
 *
 *  @copydetails SdssShapeControl
 */
class SdssShapeControl {
public:
    LSST_CONTROL_FIELD(background, double, "Additional value to add to background");
    LSST_CONTROL_FIELD(maxIter, int, "Maximum number of iterations");
    LSST_CONTROL_FIELD(maxShift, double, "Maximum centroid shift, limited to 2-10");
    LSST_CONTROL_FIELD(tol1, float, "Convergence tolerance for e1,e2");
    LSST_CONTROL_FIELD(tol2, float, "Convergence tolerance for FWHM");
    LSST_CONTROL_FIELD(doMeasurePsf, bool, "Whether to also compute the shape of the PSF model");

    /// @copydoc SdssShapeControl::SdssShapeControl
    SdssShapeControl() : background(0.0), maxIter(100), maxShift(), tol1(1E-5), tol2(1E-4),
                         doMeasurePsf(true) {}
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
     *  @param[in,out] schema        Schema to add fields to.
     *  @param[in]     name          Name prefix for all fields; "_xx", "_yy", etc. will be appended to this
     *                               to form the full field names.
     *  @param[in]     numFlags      Integer to accommodate not setting the Psf shape fields when
     *                               doMeasurePsf is false.
     *  @param[in]     doMeasurePsf  Boolean indicating whether or not the Psf is being measured (as
     *                               set in the SdssShapeControl class).
     */
    static SdssShapeResultKey addFields(
        afw::table::Schema & schema,
        std::string const & name,
        bool doMeasurePsf
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

    /// Get an SdssShapeResult from the given record
    virtual SdssShapeResult get(afw::table::BaseRecord const & record) const;

    /// Set an SdssShapeResult in the given record
    virtual void set(afw::table::BaseRecord & record, SdssShapeResult const & value) const;

    /// Get a Quadrupole for the Psf from the given record
    virtual afw::geom::ellipses::Quadrupole getPsfShape(afw::table::BaseRecord const & record) const;

    /// Set a Quadrupole for the Psf at the position of the given record
    virtual void setPsfShape(afw::table::BaseRecord & record,
                             afw::geom::ellipses::Quadrupole const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Keys
    bool operator==(SdssShapeResultKey const & other) const;
    bool operator!=(SdssShapeResultKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if the key is valid.
    bool isValid() const;

    FlagHandler const & getFlagHandler() const { return _flagHandler; }

private:
    bool _includePsf;
    ShapeResultKey _shapeResult;
    CentroidResultKey _centroidResult;
    FluxResultKey _fluxResult;
    afw::table::QuadrupoleKey _psfShapeResult;
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
 *
 *  See Bernstein & Jarvis, 2002, for more information on this type of algorithm.   Note that the
 *  code here makes no attempt to correct for the PSF; for PSF corrected ellipticities using
 *  weighted moments please use the shapeHSM package.
 */
class SdssShapeAlgorithm : public SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    static FlagDefinitionList const & getFlagDefinitions();
    static unsigned int const N_FLAGS = 6;
    static FlagDefinition const FAILURE;
    static FlagDefinition const UNWEIGHTED_BAD;
    static FlagDefinition const UNWEIGHTED;
    static FlagDefinition const SHIFT;
    static FlagDefinition const MAXITER;
    static FlagDefinition const PSF_SHAPE_BAD;

    typedef SdssShapeControl Control;
    typedef SdssShapeResult Result;
    typedef SdssShapeResultKey ResultKey;

    // NOTE: In order to accommodate the optional setting of additional fields when running with
    //       doMeasurePsf = true (do set extra fields) or false (do NOT set extra fields), all of
    //       the code in SdssShape assumes that PSF_SHAPE_BAD is the last entry in the enum list.
    //       If new flags are added, be sure to add them above the PSF_SHAPE_BAD entry.

    SdssShapeAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    /**
     *  Compute the adaptive Gaussian-weighted moments of an image.
     *
     *  @param[in] image     An Image or MaskedImage instance with int, float, or double pixels.  This
     *                       need not be a small postage stamp (the pixel region actually used in the
     *                       fit will be a subset of this image determined automatically).
     *  @param[in] position  Center position of the object to be measured, in the image's PARENT coordinates.
     *  @param[in] negative  Boolean, specify if the source is in negative flux space
     *  @param[in] ctrl      Control object specifying the details of how the object is to be measured.
     */
    template <typename ImageT>
    static Result computeAdaptiveMoments(
        ImageT const & image,
        afw::geom::Point2D const & position,
        bool negative=false,
        Control const & ctrl=Control()
    );

    /**
     *  Compute the flux within a fixed Gaussian aperture.
     *
     *  @param[in] image     An Image or MaskedImage instance with int, float, or double pixels.  This
     *                       need not be a small postage stamp (the pixel region actually used in the
     *                       fit will be a subset of this image determined automatically).
     *  @param[in] shape     Ellipse object specifying the 1-sigma contour of the Gaussian.
     *  @param[in] position  Center position of the object to be measured, in the image's PARENT coordinates.
     */
    template <typename ImageT>
    static FluxResult computeFixedMomentsFlux(
        ImageT const & image,
        afw::geom::ellipses::Quadrupole const & shape,
        afw::geom::Point2D const & position
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
    ErrElement flux_xx_Cov; ///< flux, xx term in the uncertainty covariance matrix
    ErrElement flux_yy_Cov; ///< flux, yy term in the uncertainty covariance matrix
    ErrElement flux_xy_Cov; ///< flux, xy term in the uncertainty covariance matrix

#ifndef SWIG
    std::bitset<SdssShapeAlgorithm::N_FLAGS> flags; ///< Status flags (see SdssShapeAlgorithm).
#endif

    /// Flag getter for Swig, which doesn't understand std::bitset
    bool getFlag(unsigned int index) { return flags[index]; }

    bool getFlag(std::string name) {
       return flags[SdssShapeAlgorithm::getFlagDefinitions().getDefinition(name).number];
    }

    SdssShapeResult(); ///< Constructor; initializes everything to NaN

};

/**
 *  Transformation for SdssShape measurements.
 *
 *  SdssShape measures not just shape but also flux and centroid. This
 *  transform operates on the first directly, and delegates to the Flux and
 *  Centroid transforms for the other two.
 */
class SdssShapeTransform : public BaseTransform {
public:
    typedef SdssShapeControl Control;

    SdssShapeTransform(Control const & ctrl, std::string const & name, afw::table::SchemaMapper & mapper);

    /*
     * @brief Perform transformation from inputCatalog to outputCatalog.
     *
     * @param[in]     inputCatalog   Source of data to be transformed
     * @param[in,out] outputCatalog  Container for transformed results
     * @param[in]     wcs            World coordinate system under which transformation will take place
     * @param[in]     calib          Photometric calibration under which transformation will take place
     * @throws        LengthError    Catalog sizes do not match
     */
    virtual void operator()(afw::table::SourceCatalog const & inputCatalog,
                            afw::table::BaseCatalog & outputCatalog,
                            afw::image::Wcs const & wcs,
                            afw::image::Calib const & calib) const;
private:
    FluxTransform _fluxTransform;
    CentroidTransform _centroidTransform;
    ShapeResultKey _outShapeKey;
    afw::table::QuadrupoleKey _outPsfShapeKey;
    bool _transformPsf;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssShape_h_INCLUDED
