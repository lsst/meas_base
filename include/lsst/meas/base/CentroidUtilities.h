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

#ifndef LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED
#define LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED

#include "lsst/meas/base/constants.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/table/aggregates.h"
#include "lsst/meas/base/Transform.h"

namespace lsst {
namespace meas {
namespace base {

/**
 *  @brief A reusable struct for centroid measurements.
 */
struct CentroidResult {
    CentroidElement x;   ///< x (column) coordinate of the measured position
    CentroidElement y;   ///< y (row) coordinate of the measured position
    ErrElement xErr;     ///< standard deviation of x
    ErrElement yErr;     ///< standard deviation of y
    ErrElement x_y_Cov;  ///< x,y term in the uncertainty convariance matrix

    /// Constructor; initializes everything to NaN.
    CentroidResult();

    /// Constructor; initializes everything from values.
    explicit CentroidResult(CentroidElement x_, CentroidElement y_, CentroidCov const& matrix)
            : x(x_),
              y(y_),
              xErr(std::sqrt(matrix(0, 0))),
              yErr(std::sqrt(matrix(1, 1))),
              x_y_Cov(matrix(0, 1)) {}

    /// Constructor; initializes everything from values.
    explicit CentroidResult(CentroidElement x_, CentroidElement y_, ErrElement xErr_, ErrElement yErr_)
            : x(x_), y(y_), xErr(xErr_), yErr(yErr_), x_y_Cov(0.0) {}

    /// Return a Point object containing the measured x and y
    Centroid const getCentroid() const;

    /// Set the struct fields from the given Point object.
    void setCentroid(Centroid const& centroid);

    /// Return the 2D point type corresponding to this result
    geom::Point<CentroidElement> getPoint() { return geom::Point<CentroidElement>(x, y); }

    /// Return the 2x2 symmetric covariance matrix, with rows and columns ordered (x, y)
    CentroidCov const getCentroidErr() const;

    /// Set the struct uncertainty fields from the given matrix, with rows and columns ordered (x, y)
    void setCentroidErr(CentroidCov const& matrix);

    /// Set the struct uncertainty fields from the sigma values
    void setCentroidErr(ErrElement _xErr, ErrElement _yErr);
};

std::ostream& operator<<(std::ostream& os, CentroidResult const& result);

/**
 *  @brief A FunctorKey for CentroidResult
 *
 *  This class makes it easy to copy centroids and their uncertainties to and from records, and provides
 *  a method to add the appropriate fields to a Schema.
 */
class CentroidResultKey : public afw::table::FunctorKey<CentroidResult> {
public:
    /**
     *  @brief Add the appropriate fields to a Schema, and return a CentroidResultKey that manages them
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_x", "_y", etc. will be appended to this
     *                         to form the full field names.
     *  @param[in]     doc     String used as the documentation for the x and y fields.
     *  @param[in] uncertainty Enum indicating which (if any) uncertainty values will be saved.
     *
     *  The unit for all fields will be pixels (pixels^2 for covariances).
     */
    static CentroidResultKey addFields(afw::table::Schema& schema, std::string const& name,
                                       std::string const& doc, UncertaintyEnum uncertainty);

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    CentroidResultKey() : _centroid(), _centroidErr() {}

    /// Construct from a pair of Keys
    CentroidResultKey(afw::table::PointKey<CentroidElement> const& centroid,
                      afw::table::CovarianceMatrixKey<ErrElement, 2> const& centroidErr)
            : _centroid(centroid), _centroidErr(centroidErr) {}

    /**
     *  @brief Construct from a subschema, assuming _x, _y, etc. subfields
     *
     *  If a schema has "a_x", "a_y", etc. fields, this constructor allows you to construct
     *  a CentroidResultKey via:
     *  @code
     *  CentroidResultKey k(schema["a"]);
     *  @endcode
     */
    CentroidResultKey(afw::table::SubSchema const& s);

    /// Get a CentroidResult from the given record
    virtual CentroidResult get(afw::table::BaseRecord const& record) const;

    /// Set a CentroidResult in the given record
    virtual void set(afw::table::BaseRecord& record, CentroidResult const& value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Keys
    bool operator==(CentroidResultKey const& other) const {
        return _centroid == other._centroid && _centroidErr == other._centroidErr;
    }
    bool operator!=(CentroidResultKey const& other) const { return !(*this == other); }
    //@}

    /// Return True if the centroid key is valid.
    bool isValid() const { return _centroid.isValid() && _centroidErr.isValid(); }

    /// Return a FunctorKey to just the centroid value
    afw::table::PointKey<CentroidElement> getCentroid() const { return _centroid; }

    /// Return a FunctorKey to just the uncertainty matrix
    afw::table::CovarianceMatrixKey<ErrElement, 2> getCentroidErr() const { return _centroidErr; }

    /// Return a Key for the x coordinate
    afw::table::Key<CentroidElement> getX() const { return _centroid.getX(); }

    /// Return a Key for the y coordinate
    afw::table::Key<CentroidElement> getY() const { return _centroid.getY(); }

private:
    afw::table::PointKey<CentroidElement> _centroid;
    afw::table::CovarianceMatrixKey<ErrElement, 2> _centroidErr;
};

/**
 *  Base for centroid measurement transformations.
 *
 *  Provides a basic transform from centroid plus associated uncertainty to
 *  celestial position with uncertainty. The basic "flag" attribute for the
 *  measurement algorithm is propagated to the output.
 *
 *  Subclasses should define a constructor which take a Control argument
 *  corresponding to the measurement algorithm being transformed and ensure
 *  that any other necessary flags are propagated.
 */
class CentroidTransform : public BaseTransform {
public:
    CentroidTransform(std::string const& name, afw::table::SchemaMapper& mapper);

    /*
     * @brief Perform transformation from inputCatalog to outputCatalog.
     *
     * @param[in]     inputCatalog   Source of data to be transformed
     * @param[in,out] outputCatalog  Container for transformed results
     * @param[in]     wcs            World coordinate system under which transformation will take place
     * @param[in]     photoCalib     Photometric calibration under which transformation will take place
     * @throws        LengthError    Catalog sizes do not match
     */
    virtual void operator()(afw::table::SourceCatalog const& inputCatalog,
                            afw::table::BaseCatalog& outputCatalog, afw::geom::SkyWcs const& wcs,
                            afw::image::PhotoCalib const& photoCalib) const;

private:
    afw::table::CoordKey _coordKey;
    afw::table::CovarianceMatrixKey<ErrElement, 2> _coordErrKey;
};

class CentroidChecker {
public:
    /**
     *  Check source record produced by a centroid algorithm called "name".
     *
     *  If the centroid is accompanied by uncertainties (the xErr and yErr
     *  fields), these should not be NaN. If they are, the algorithmName +
     *  "_flag_badError" flag and the general failure flag for the algorithm
     *  are both set.
     *
     *  If the centroid set by the algorithm lies outside the footprint attached
     *  to the record, or the centroid is more than "dist" pixels from the
     *  footprint peak:
     *
     *     (1) the general failure flag for the algorithm is set
     *     (2) algorithmName + "_flag_resetToPeak" flag is set
     *     (3) the position of the centroid is changed to that of the footprint Peak
     *
     *  @param[in,out] schema        Schema to which the flag_resetToPeak is to be added
     *  @param[in]  name             The name of the algorithm we will be checking
     *  @param[in]  doFootprintCheck Check if centroid is within footprint
     *  @param[in]  maxDistFromPeak  If >0; maximum distance in pixels between the footprint
     *                                peak and centroid allowed before resetToPeak flag is set.
     */
    CentroidChecker(afw::table::Schema& schema, std::string const& name, bool inside = true,
                    double maxDistFromPeak = -1.0);

    /**
     *  Set the centroid to the first footprint if the centroid is either more than _dist
     *  pixels from the footprint center, or if it is outside the footprint. Set appropriate
     *  flags to indicate any changes to the centroid, and to indicate if uncertainties
     *  are set to invalid ("NaN") values.
     */
    bool operator()(afw::table::SourceRecord& record) const;

private:
    bool _doFootprintCheck;
    double _maxDistFromPeak;
    afw::table::Key<afw::table::Flag> _resetKey;
    afw::table::Key<afw::table::Flag> _failureKey;
    afw::table::Key<afw::table::Flag> _badErrorKey;
    afw::table::Key<CentroidElement> _xKey;
    afw::table::Key<CentroidElement> _yKey;
    afw::table::Key<ErrElement> _xErrKey;
    afw::table::Key<ErrElement> _yErrKey;
};
}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED
