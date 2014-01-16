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

#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

BaseAlgorithmMapper::BaseAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix
) :
    _flag(
        schema.addField(
            afw::table::Field<afw::table::Flag>(
                prefix + "_flag",
                "general failure flag for " + prefix + " measurement"
            ),
            true // replace existing fields if present
        )
    )
{}

void BaseAlgorithmMapper::fail(afw::table::BaseRecord & record) {
    record.set(_flag, true);
}

FluxAlgorithmMapper::FluxAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    ResultMapperUncertaintyEnum uncertainty
) : BaseAlgorithmMapper(schema, prefix) {
    _flux = schema.addField(
        afw::table::Field<Flux>(prefix + "_flux", "measured flux", "dn"),
        true
    );
    if (uncertainty != NO_UNCERTAINTY) {
        _fluxSigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_fluxSigma", "1-sigma uncertainty on " + prefix + "_flux", "dn"
            ),
            true
         );
    }
}

void FluxAlgorithmMapper::apply(afw::table::BaseRecord & record, FluxAlgorithmResult const & result) {
    record.set(_flux, result.flux);
    if (_fluxSigma.isValid()) {
        record.set(_fluxSigma, result.fluxSigma);
    }
    record.set(_flag, false);
}

CentroidAlgorithmMapper::CentroidAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    ResultMapperUncertaintyEnum uncertainty
) : BaseAlgorithmMapper(schema, prefix) {
    _x = schema.addField(
        afw::table::Field<CentroidElement>(
            prefix + "_x",
            "x coordinate of position",
            "pixels"
        ),
        true // replace existing fields if present
    );
    _y = schema.addField(
        afw::table::Field<CentroidElement>(
            prefix + "_y",
            "y coordinate of position",
            "pixels"
        ),
        true // replace existing fields if present
    );
    if (uncertainty != NO_UNCERTAINTY) {
        _xSigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_xSigma", "1-sigma uncertainty on " + prefix + "_x", "pixels"
            ),
            true
         );
        _ySigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_ySigma", "1-sigma uncertainty on " + prefix + "_y", "pixels"
            ),
            true
        );
        if (uncertainty == FULL_COVARIANCE) {
            _x_y_Cov = schema.addField(
                afw::table::Field<ErrElement>(
                    prefix + "_x_y_Cov",
                    "uncertainty covariance between " + prefix + "_x and " + prefix + "_y",
                    "pixels^2"
                ),
                true
            );
        }
    }
}

void CentroidAlgorithmMapper::apply(afw::table::BaseRecord & record, CentroidAlgorithmResult const & result) {
    record.set(_x, result.x);
    record.set(_y, result.y);
    if (_xSigma.isValid()) {
        assert(_ySigma.isValid());
        record.set(_xSigma, result.xSigma);
        record.set(_ySigma, result.ySigma);
        if (_x_y_Cov.isValid()) {
            record.set(_x_y_Cov, result.x_y_Cov);
        }
    }
    record.set(_flag, false);
}

ShapeAlgorithmMapper::ShapeAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    ResultMapperUncertaintyEnum uncertainty
) : BaseAlgorithmMapper(schema, prefix) {
    _xx = schema.addField(
        afw::table::Field<ShapeElement>(
            prefix + "_xx",
            "x-x second moment of ellipse",
            "pixels^2"
        ),
        true // replace existing fields if present
    );
    _yy = schema.addField(
        afw::table::Field<ShapeElement>(
            prefix + "_yy",
            "y-y second moment of ellipse",
            "pixels^2"
        ),
        true // replace existing fields if present
    );
    _xy = schema.addField(
        afw::table::Field<ShapeElement>(
            prefix + "_xy",
            "x-y second moment of ellipse",
            "pixels^2"
        ),
        true // replace existing fields if present
    );
    if (uncertainty != NO_UNCERTAINTY) {
        _xxSigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_xxSigma", "1-sigma uncertainty on " + prefix + "_xx", "pixels^2"
            ),
            true
         );
        _yySigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_yySigma", "1-sigma uncertainty on " + prefix + "_yy", "pixels^2"
            ),
            true
        );
        _xySigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_xySigma", "1-sigma uncertainty on " + prefix + "_xy", "pixels^2"
            ),
            true
        );
        if (uncertainty == FULL_COVARIANCE) {
            _xx_yy_Cov = schema.addField(
                afw::table::Field<ErrElement>(
                    prefix + "_xx_yy_Cov",
                    "uncertainty covariance between " + prefix + "_xx and " + prefix + "_yy",
                    "pixels^4"
                ),
                true
            );
            _xx_xy_Cov = schema.addField(
                afw::table::Field<ErrElement>(
                    prefix + "_xx_xy_Cov",
                    "uncertainty covariance between " + prefix + "_xx and " + prefix + "_xy",
                    "pixels^4"
                ),
                true
            );
            _yy_xy_Cov = schema.addField(
                afw::table::Field<ErrElement>(
                    prefix + "_yy_xy_Cov",
                    "uncertainty covariance between " + prefix + "_yy and " + prefix + "_xy",
                    "pixels^4"
                ),
                true
            );
        }
    }
}

void ShapeAlgorithmMapper::apply(afw::table::BaseRecord & record, ShapeAlgorithmResult const & result) {
    record.set(_xx, result.xx);
    record.set(_yy, result.yy);
    record.set(_xy, result.xy);
    if (_xxSigma.isValid()) {
        assert(_yySigma.isValid());
        assert(_xySigma.isValid());
        record.set(_xxSigma, result.xxSigma);
        record.set(_yySigma, result.yySigma);
        record.set(_xySigma, result.xySigma);
        if (_xx_yy_Cov.isValid()) {
            assert(_xx_xy_Cov.isValid());
            assert(_yy_xy_Cov.isValid());
            record.set(_xx_yy_Cov, result.xx_yy_Cov);
            record.set(_xx_xy_Cov, result.xx_xy_Cov);
            record.set(_yy_xy_Cov, result.yy_xy_Cov);
        }
    }
    record.set(_flag, false);
}

}}} // namespace lsst::meas::base
