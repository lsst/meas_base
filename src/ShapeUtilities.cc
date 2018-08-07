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

#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/afw/table/BaseRecord.h"

namespace lsst {
namespace meas {
namespace base {

ShapeResult::ShapeResult()
        : xx(std::numeric_limits<ShapeElement>::quiet_NaN()),
          yy(std::numeric_limits<ShapeElement>::quiet_NaN()),
          xy(std::numeric_limits<ShapeElement>::quiet_NaN()),
          xxErr(std::numeric_limits<ErrElement>::quiet_NaN()),
          yyErr(std::numeric_limits<ErrElement>::quiet_NaN()),
          xyErr(std::numeric_limits<ErrElement>::quiet_NaN()),
          xx_yy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
          xx_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
          yy_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()) {}

Shape const ShapeResult::getShape() const { return Shape(xx, yy, xy); }

void ShapeResult::setShape(Shape const &shape) {
    xx = shape.getIxx();
    yy = shape.getIyy();
    xy = shape.getIxy();
}

ShapeCov const ShapeResult::getShapeErr() const {
    ShapeCov m;
    m << xxErr * xxErr, xx_yy_Cov, xx_xy_Cov, xx_yy_Cov, yyErr * yyErr, yy_xy_Cov, xx_xy_Cov,
            yy_xy_Cov, xyErr * xyErr;
    return m;
}

void ShapeResult::setShapeErr(ShapeCov const &matrix) {
    xxErr = std::sqrt(matrix(0, 0));
    yyErr = std::sqrt(matrix(1, 1));
    xyErr = std::sqrt(matrix(2, 2));
    xx_yy_Cov = matrix(0, 1);
    xx_xy_Cov = matrix(0, 2);
    yy_xy_Cov = matrix(1, 2);
}

void ShapeResult::setShapeErr(ErrElement _xxErr, ErrElement _yyErr, ErrElement _xyErr) {
    xxErr = _xxErr;
    yyErr = _yyErr;
    xyErr = _xyErr;
    xx_yy_Cov = 0.0;
    xx_xy_Cov = 0.0;
    yy_xy_Cov = 0.0;
}

ShapeResultKey ShapeResultKey::addFields(afw::table::Schema &schema, std::string const &name,
                                         std::string const &doc, UncertaintyEnum uncertainty,
                                         afw::table::CoordinateType coordType) {
    ShapeResultKey r;
    r._shape = afw::table::QuadrupoleKey::addFields(schema, name, doc, coordType);
    if (uncertainty != NO_UNCERTAINTY) {
        std::vector<afw::table::Key<ErrElement> > sigma(3);
        std::vector<afw::table::Key<ErrElement> > cov;
        sigma[0] = schema.addField<ErrElement>(
                schema.join(name, "xxErr"), "Standard deviation of xx moment",
                coordType == afw::table::CoordinateType::PIXEL ? "pixel^2" : "rad^2");
        sigma[1] = schema.addField<ErrElement>(
                schema.join(name, "yyErr"), "Standard deviation of yy moment",
                coordType == afw::table::CoordinateType::PIXEL ? "pixel^2" : "rad^2");
        sigma[2] = schema.addField<ErrElement>(
                schema.join(name, "xyErr"), "Standard deviation of xy moment",
                coordType == afw::table::CoordinateType::PIXEL ? "pixel^2" : "rad^2");
        if (uncertainty == FULL_COVARIANCE) {
            cov.push_back(schema.addField<ErrElement>(
                    schema.join(name, "xx_yy_Cov"), "uncertainty covariance in xx and yy",
                    coordType == afw::table::CoordinateType::PIXEL ? "pixel^4" : "rad^4"));
            cov.push_back(schema.addField<ErrElement>(
                    schema.join(name, "xx_xy_Cov"), "uncertainty covariance in xx and xy",
                    coordType == afw::table::CoordinateType::PIXEL ? "pixel^4" : "rad^4"));
            cov.push_back(schema.addField<ErrElement>(
                    schema.join(name, "yy_xy_Cov"), "uncertainty covariance in yy and xy",
                    coordType == afw::table::CoordinateType::PIXEL ? "pixel^4" : "rad^4"));
        }
        r._shapeErr = afw::table::CovarianceMatrixKey<ErrElement, 3>(sigma, cov);
    }
    return r;
}

namespace {

std::vector<std::string> getNameVector() {
    std::vector<std::string> v;
    v.push_back("xx");
    v.push_back("yy");
    v.push_back("xy");
    return v;
}

}  // namespace

ShapeResultKey::ShapeResultKey(afw::table::SubSchema const &s) : _shape(s) {
    static std::vector<std::string> names = getNameVector();  // C++11 TODO: just use initializer list
    try {
        _shapeErr = afw::table::CovarianceMatrixKey<ErrElement, 3>(s, names);
    } catch (pex::exceptions::NotFoundError &) {
    }
}

ShapeResult ShapeResultKey::get(afw::table::BaseRecord const &record) const {
    ShapeResult r;
    r.setShape(record.get(_shape));
    if (_shapeErr.isValid()) {
        r.setShapeErr(record.get(_shapeErr));
    }
    return r;
}

void ShapeResultKey::set(afw::table::BaseRecord &record, ShapeResult const &value) const {
    record.set(_shape, value.getShape());
    if (_shapeErr.isValid()) {
        record.set(_shapeErr, value.getShapeErr());
    }
}

ShapeTrMatrix makeShapeTransformMatrix(geom::LinearTransform const &xform) {
    typedef geom::LinearTransform LT;
    Eigen::Matrix<ShapeElement, 3, 3, Eigen::DontAlign> m;
    m << xform[LT::XX] * xform[LT::XX], xform[LT::XY] * xform[LT::XY], 2 * xform[LT::XX] * xform[LT::XY],
            xform[LT::YX] * xform[LT::YX], xform[LT::YY] * xform[LT::YY], 2 * xform[LT::YX] * xform[LT::YY],
            xform[LT::XX] * xform[LT::YX], xform[LT::XY] * xform[LT::YY],
            xform[LT::XX] * xform[LT::YY] + xform[LT::XY] * xform[LT::YX];
    return m;
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
