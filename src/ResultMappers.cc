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

template <std::size_t N>
FlagsResultMapper<N>::FlagsResultMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    boost::array<FlagDef,N> const & flagDefs
) {
    _flags[0] = schema.addField(
        afw::table::Field<afw::table::Flag>(
            prefix + "_flag",
            "general failure flag for " + prefix + " measurement"
        ),
        true // replace existing fields if present
    );
    for (std::size_t i = 0; i < N; ++i) {
        _flags[i+1] = schema.addField(
            afw::table::Field<afw::table::Flag>(flagDefs[i].name, flagDefs[i].doc),
            true // replace existing fields if present
        );
    }
}

template <std::size_t N>
void FlagsResultMapper<N>::fail(afw::table::BaseRecord & record, MeasurementError const & error) const {
    assert(error.getFlagBit() < N);
    record.set(_flags[0], true);
    record.set(_flags[error.getFlagBit() + 1], true);
}

template <std::size_t N>
void FlagsResultMapper<N>::apply(afw::table::BaseRecord & record, FlagsResult<N> const & result) const {
    for (std::size_t i = 0; i < N; ++i) {
        record.set(_flags[i+1], result.flags[i]);
    }
    record.set(_flags[0], false);
}

FluxResultMapper::FluxResultMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    ResultMapperUncertaintyEnum uncertainty
) {
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

void FluxResultMapper::apply(afw::table::BaseRecord & record, FluxResult const & result) const {
    record.set(_flux, result.flux);
    if (_fluxSigma.isValid()) {
        record.set(_fluxSigma, result.fluxSigma);
    }
}

CentroidResultMapper::CentroidResultMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    ResultMapperUncertaintyEnum uncertainty
) {
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

void CentroidResultMapper::apply(afw::table::BaseRecord & record, CentroidResult const & result) const {
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
}

ShapeResultMapper::ShapeResultMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    ResultMapperUncertaintyEnum uncertainty
) {
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

void ShapeResultMapper::apply(afw::table::BaseRecord & record, ShapeResult const & result) const {
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
}

template class FlagsResultMapper<0>;
template class FlagsResultMapper<1>;
template class FlagsResultMapper<2>;
template class FlagsResultMapper<3>;
template class FlagsResultMapper<4>;
template class FlagsResultMapper<5>;
template class FlagsResultMapper<6>;
template class FlagsResultMapper<7>;
template class FlagsResultMapper<8>;

}}} // namespace lsst::meas::base
