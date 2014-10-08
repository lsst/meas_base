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

#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/afw/table/BaseRecord.h"

namespace lsst { namespace meas { namespace base {

ShapeResult::ShapeResult() :
    xx(std::numeric_limits<ShapeElement>::quiet_NaN()),
    yy(std::numeric_limits<ShapeElement>::quiet_NaN()),
    xy(std::numeric_limits<ShapeElement>::quiet_NaN()),
    xxSigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    yySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    xySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    xx_yy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    xx_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    yy_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}

Shape const ShapeResult::getShape() const { return Shape(xx, yy, xy); }

void ShapeResult::setShape(Shape const & shape) {
    xx = shape.getIxx();
    yy = shape.getIyy();
    xy = shape.getIxy();
}

ShapeCov const ShapeResult::getShapeErr() const {
    ShapeCov m;
    m <<
        xxSigma*xxSigma, xx_yy_Cov, xx_xy_Cov,
        xx_yy_Cov, yySigma*yySigma, yy_xy_Cov,
        xx_xy_Cov, yy_xy_Cov, xySigma*xySigma;
    return m;
}

void ShapeResult::setShapeErr(ShapeCov const & matrix) {
    xxSigma = std::sqrt(matrix(0, 0));
    yySigma = std::sqrt(matrix(1, 1));
    xySigma = std::sqrt(matrix(2, 2));
    xx_yy_Cov = matrix(0, 1);
    xx_xy_Cov = matrix(0, 2);
    yy_xy_Cov = matrix(1, 2);
}

ShapeResultKey ShapeResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc,
    UncertaintyEnum uncertainty
) {
    ShapeResultKey r;
    r._shape = afw::table::QuadrupoleKey::addFields(
        schema,
        name,
        doc,
        "pixels^2"
    );
    if (uncertainty != NO_UNCERTAINTY) {
        std::vector< afw::table::Key<ErrElement> > sigma(3);
        std::vector< afw::table::Key<ErrElement> > cov;
        sigma[0] = schema.addField<ErrElement>(
            schema.join(name, "xxSigma"), "1-sigma uncertainty on xx moment", "pixels^2"
        );
        sigma[1] = schema.addField<ErrElement>(
            schema.join(name, "yySigma"), "1-sigma uncertainty on yy moment", "pixels^2"
        );
        sigma[2] = schema.addField<ErrElement>(
            schema.join(name, "xySigma"), "1-sigma uncertainty on xy moment", "pixels^2"
        );
        if (uncertainty == FULL_COVARIANCE) {
            cov.push_back(
                schema.addField<ErrElement>(
                    schema.join(name, "xx_yy_Cov"), "uncertainty covariance in xx and yy", "pixels^4"
                )
            );
            cov.push_back(
                schema.addField<ErrElement>(
                    schema.join(name, "xx_xy_Cov"), "uncertainty covariance in xx and xy", "pixels^4"
                )
            );
            cov.push_back(
                schema.addField<ErrElement>(
                    schema.join(name, "yy_xy_Cov"), "uncertainty covariance in yy and xy", "pixels^4"
                )
            );
        }
        r._shapeErr = afw::table::CovarianceMatrixKey<ErrElement,3>(sigma, cov);
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

} // anonymous

ShapeResultKey::ShapeResultKey(afw::table::SubSchema const & s) :
    _shape(s)
{
    static std::vector<std::string> names = getNameVector(); // C++11 TODO: just use initializer list
    try {
        _shapeErr = afw::table::CovarianceMatrixKey<ErrElement,3>(s, names);
    } catch (pex::exceptions::NotFoundError &) {}
}

ShapeResult ShapeResultKey::get(afw::table::BaseRecord const & record) const {
    ShapeResult r;
    r.setShape(record.get(_shape));
    if (_shapeErr.isValid()) {
        r.setShapeErr(record.get(_shapeErr));
    }
    return r;
}

void ShapeResultKey::set(afw::table::BaseRecord & record, ShapeResult const & value) const {
    record.set(_shape, value.getShape());
    if (_shapeErr.isValid()) {
        record.set(_shapeErr, value.getShapeErr());
    }
}

}}} // lsst::meas::base
