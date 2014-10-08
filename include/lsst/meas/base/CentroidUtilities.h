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
#include "lsst/afw/table/aggregates.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A reusable struct for centroid measurements.
 */
struct CentroidResult {
    CentroidElement x; ///< x (column) coordinate of the measured position
    CentroidElement y; ///< y (row) coordinate of the measured position
    ErrElement xSigma; ///< 1-Sigma uncertainty on x (sqrt of variance)
    ErrElement ySigma; ///< 1-Sigma uncertainty on y (sqrt of variance)
    ErrElement x_y_Cov; ///< x,y term in the uncertainty convariance matrix

    /// Constructor; initializes everything to NaN.
    CentroidResult();

    /// Return a Point object containing the measured x and y
    Centroid const getCentroid() const;

    /// Set the struct fields from the given Point object.
    void setCentroid(Centroid const & centroid);

    /// Return the 2x2 symmetric covariance matrix, with rows and columns ordered (x, y)
    CentroidCov const getCentroidErr() const;

    /// Set the struct uncertainty fields from the given matrix, with rows and columns ordered (x, y)
    void setCentroidErr(CentroidCov const & matrix);

};

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
    static CentroidResultKey addFields(
        afw::table::Schema & schema,
        std::string const & name,
        std::string const & doc,
        UncertaintyEnum uncertainty
    );

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    CentroidResultKey() : _centroid(), _centroidErr() {}

    /// Construct from a pair of Keys
    CentroidResultKey(
        afw::table::PointKey<CentroidElement> const & centroid,
        afw::table::CovarianceMatrixKey<ErrElement,2> const & centroidErr
    ) :
        _centroid(centroid), _centroidErr(centroidErr)
    {}

    /**
     *  @brief Construct from a subschema, assuming _x, _y, etc. subfields
     *
     *  If a schema has "a_x", "a_y", etc. fields, this constructor allows you to construct
     *  a CentroidResultKey via:
     *  @code
     *  CentroidResultKey k(schema["a"]);
     *  @endcode
     */
    CentroidResultKey(afw::table::SubSchema const & s);

    /// Get a CentroidResult from the given record
    virtual CentroidResult get(afw::table::BaseRecord const & record) const;

    /// Set a CentroidResult in the given record
    virtual void set(afw::table::BaseRecord & record, CentroidResult const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Keys
    bool operator==(CentroidResultKey const & other) const {
        return _centroid == other._centroid && _centroidErr == other._centroidErr;
    }
    bool operator!=(CentroidResultKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if the centroid key is valid.
    bool isValid() const { return _centroid.isValid() && _centroidErr.isValid(); }

    /// Return a FunctorKey to just the centroid value
    afw::table::PointKey<CentroidElement> getCentroid() const { return _centroid; }

    /// Return a FunctorKey to just the uncertainty matrix
    afw::table::CovarianceMatrixKey<ErrElement,2> getCentroidErr() const { return _centroidErr; }

    /// Return a Key for the x coordinate
    afw::table::Key<CentroidElement> getX() const { return _centroid.getX(); }

    /// Return a Key for the y coordinate
    afw::table::Key<CentroidElement> getY() const { return _centroid.getY(); }

private:
    afw::table::PointKey<CentroidElement> _centroid;
    afw::table::CovarianceMatrixKey<ErrElement,2> _centroidErr;
};

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_CentroidUtilities_h_INCLUDED
