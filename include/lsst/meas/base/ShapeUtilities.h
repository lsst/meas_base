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

#ifndef LSST_MEAS_BASE_ShapeUtilities_h_INCLUDED
#define LSST_MEAS_BASE_ShapeUtilities_h_INCLUDED

#include "lsst/meas/base/constants.h"
#include "lsst/afw/table/aggregates.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A reusable struct for moments-based shape measurements.
 *
 *  Shape measurements and their errors should always be in pixels coordinates.  This struct should generally
 *  be preferred over a custom struct with other ellipse parametrizations unless the measurement takes place
 *  in another parametrization and a transformation to this one would result in a loss of information or
 *  obfuscate the results of the measurement (i.e. use this one unless you have a good reason not to).
 */
struct ShapeResult {
    ShapeElement xx; // image or model second moment for x^2
    ShapeElement yy; // image or model second moment for y^2
    ShapeElement xy; // image or model second moment for xy^2
    ErrElement xxSigma; ///< 1-Sigma uncertainty on xx (sqrt of variance)
    ErrElement yySigma; ///< 1-Sigma uncertainty on yy (sqrt of variance)
    ErrElement xySigma; ///< 1-Sigma uncertainty on xy (sqrt of variance)
    ErrElement xx_yy_Cov; ///< xx,yy term in the uncertainty convariance matrix
    ErrElement xx_xy_Cov; ///< xx,xy term in the uncertainty convariance matrix
    ErrElement yy_xy_Cov; ///< yy,xy term in the uncertainty convariance matrix

    /// Constructor; initializes everything to NaN.
    ShapeResult();

    /// Constructor; initializes everything from values.
    explicit ShapeResult(ShapeElement _xx, ShapeElement _yy, ShapeElement _xy, ShapeCov const & matrix):
        xx(_xx),
        yy(_yy),
        xy(_xy),
        xxSigma(std::sqrt(matrix(0, 0))),
        yySigma(std::sqrt(matrix(1, 1))),
        xySigma(std::sqrt(matrix(2, 2))),
        xx_yy_Cov(matrix(0, 1)),
        xx_xy_Cov(matrix(0, 2)),
        yy_xy_Cov(matrix(1, 2))
    {}

    /// Constructor; initializes everything from values.
    explicit ShapeResult(ShapeElement _xx, ShapeElement _yy, ShapeElement _xy,
                            ErrElement _xxSigma, ErrElement _yySigma, ErrElement _xySigma) :
        xx(_xx),
        yy(_yy),
        xy(_xy),
        xxSigma(_xxSigma),
        yySigma(_yySigma),
        xySigma(_xySigma),
        xx_yy_Cov(0.0),
        xx_xy_Cov(0.0),
        yy_xy_Cov(0.0)
    {}

    /**
     *  @brief Return an afw::geom::ellipses object corresponding to xx, yy, xy.
     *
     *  This method can be used to return an average radius for the measured shape, e.g.
     *  @c getShape().getDeterminantRadius()
     */
    Shape const getShape() const;

    // get the Quadrupole corresponding to this ShapeResult
    afw::geom::ellipses::Quadrupole getQuadrupole()
    {
        return afw::geom::ellipses::Quadrupole(xx, yy, xy);
    }

    /// Set struct elements from the given Quadrupole object
    void setShape(Shape const & shape);

    /// Return the 3x3 symmetric covariance matrix, with rows and columns ordered (xx, yy, xy)
    ShapeCov const getShapeErr() const;

    /// Set the struct uncertainty elements from the given matrix, with rows and columns ordered (xx, yy, xy)
    void setShapeErr(ShapeCov const & matrix);

    /// Set the struct uncertainty elements from the given values
    void setShapeErr(ErrElement _xxSigma, ErrElement _yySigma, ErrElement _xySigma);

};

/**
 *  @brief A FunctorKey for ShapeResult
 *
 *  This class makes it easy to copy shapes and their uncertainties to and from records, and provides
 *  a method to add the appropriate fields to a Schema.
 */
class ShapeResultKey : public afw::table::FunctorKey<ShapeResult>  {
public:

    /**
     *  @brief Add the appropriate fields to a Schema, and return a ShapeResultKey that manages them
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_xx", "_yy", etc. will be appended to this
     *                         to form the full field names.
     *  @param[in]     doc     String used as the documentation for the xx, yy, xy fields.
     *  @param[in] uncertainty Enum indicating which (if any) uncertainty values will be saved.
     *
     *  The unit for all fields will be pixels^2 (pixels^4 for covariances).
     */
    static ShapeResultKey addFields(
        afw::table::Schema & schema,
        std::string const & name,
        std::string const & doc,
        UncertaintyEnum uncertainty
    );

    /// Default constructor; instance will not be usuable unless subsequently assigned to.
    ShapeResultKey() : _shape(), _shapeErr() {}

    /// Construct from a pair of Keys
    ShapeResultKey(
        afw::table::QuadrupoleKey const & shape,
        afw::table::CovarianceMatrixKey<ErrElement,3> const & shapeErr
    ) :
        _shape(shape), _shapeErr(shapeErr)
    {}

    /**
     *  @brief Construct from a subschema, assuming _xx, _yy, etc. subfields
     *
     *  If a schema has "a_xx", "a_yy", etc. fields, this constructor allows you to construct
     *  a ShapeResultKey via:
     *  @code
     *  ShapeResultKey k(schema["a"]);
     *  @endcode
     */
    ShapeResultKey(afw::table::SubSchema const & s);

    /// Get a ShapeResult from the given record
    virtual ShapeResult get(afw::table::BaseRecord const & record) const;

    /// Set a ShapeResult in the given record
    virtual void set(afw::table::BaseRecord & record, ShapeResult const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Keys
    bool operator==(ShapeResultKey const & other) const {
        return _shape == other._shape && _shapeErr == other._shapeErr;
    }
    bool operator!=(ShapeResultKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if the shape key is valid.
    bool isValid() const { return _shape.isValid() && _shapeErr.isValid(); }

    /// Return a FunctorKey to just the shape value
    afw::table::QuadrupoleKey getShape() const { return _shape; }

    /// Return a FunctorKey to just the uncertainty matrix
    afw::table::CovarianceMatrixKey<ErrElement,3> getShapeErr() const { return _shapeErr; }

    /// Return a Key for the xx moment
    afw::table::Key<ShapeElement> getIxx() const { return _shape.getIxx(); }

    /// Return a Key for the yy moment
    afw::table::Key<ShapeElement> getIyy() const { return _shape.getIyy(); }

    /// Return a Key for the xy moment
    afw::table::Key<ShapeElement> getIxy() const { return _shape.getIxy(); }

private:
    afw::table::QuadrupoleKey _shape;
    afw::table::CovarianceMatrixKey<ErrElement,3> _shapeErr;
};

/**
 *  @brief Construct a matrix suitable for transforming second moments.
 *
 *  Given an LinearTransform which maps from positions (x, y) to (a, d),
 *  returns a 3-by-3 matrix which transforms (xx, yy, xy) to (aa, dd, ad).
 *
 *  That is, given an input transform described by the matrix
 *
 *  | A11 | A12 |
 *  | A21 | A22 |
 *
 *  we return the matrix
 *
 *  | A11*A11 | A12*A12 |         2*A11*A12 |
 *  | A21*A21 | A22*A22 |         2*A21*A22 |
 *  | A11*A21 | A12*A22 | A11*A22 + A12*A21 |
 *
 *  @param[in]  xform  LinearTransform describing the coordinate mapping
 *  @returns    A 3-by-3 transformation matrix for the second order moments
 */
ShapeTrMatrix makeShapeTransformMatrix(afw::geom::LinearTransform const & xform);

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_ShapeUtilities_h_INCLUDED
