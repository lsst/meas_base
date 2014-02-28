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

#include "lsst/utils/PowFast.h"
#include "lsst/meas/base/SdssShape.h"

// Doxygen gets confused and generates warnings when trying to map the definitions here to their
// declarations, but we want to put the source code itself in the HTML docs, so we just tell it
// not to look for any documentation comments here.
/// @cond SOURCE_FILE

namespace lsst { namespace meas { namespace base {

namespace {

/*
 * The exponential function that we use, which may be only an approximation to the true value of e^x
 */
#define USE_APPROXIMATE_EXP 1
#if USE_APPROXIMATE_EXP
    lsst::utils::PowFast const& powFast = lsst::utils::getPowFast<11>();
#endif

inline float
approxExp(float x)
{
#if USE_APPROXIMATE_EXP
    return powFast.exp(x);
#else
    return std::exp(x);
#endif
}

} // anonymous

SdssShapeExtras::SdssShapeExtras() :
    xy4(std::numeric_limits<ShapeElement>::quiet_NaN()),
    xy4Sigma(std::numeric_limits<ShapeElement>::quiet_NaN()),
    flux_xx_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    flux_yy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    flux_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}

SdssShapeExtrasMapper::SdssShapeExtrasMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    SdssShapeControl const &
) :
    _xy4(
        schema.addField(
            afw::table::Field<ShapeElement>(
                // TODO: get more mathematically precise documentation on this from RHL
                prefix + "_xy4", "4th moment used in certain shear-estimation algorithms", "pixels^4"
            ), true // doReplace
        )
    ),
    _xy4Sigma(
        schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_xy4Sigma", "uncertainty on " + prefix + "_xy4", "pixels^4"
            ),
            true
        )
    ),
    _flux_xx_Cov(
        schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_flux_xx_Cov",
                "uncertainty covariance between " + prefix + "_flux and " + prefix + "_xx",
                "dn*pixels^2"
            ),
            true
        )
    ),
    _flux_yy_Cov(
        schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_flux_yy_Cov",
                "uncertainty covariance between " + prefix + "_flux and " + prefix + "_yy",
                "dn*pixels^2"
            ),
            true
        )
    ),
    _flux_xy_Cov(
        schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_flux_xy_Cov",
                "uncertainty covariance between " + prefix + "_flux and " + prefix + "_xy",
                "dn*pixels^2"
            ),
            true
        )
    )
{}

void SdssShapeExtrasMapper::apply(afw::table::BaseRecord & record, SdssShapeExtras const & result) const {
    record.set(_xy4, result.xy4);
    record.set(_xy4Sigma, result.xy4Sigma);
    record.set(_flux_xx_Cov, result.flux_xx_Cov);
    record.set(_flux_yy_Cov, result.flux_yy_Cov);
    record.set(_flux_xy_Cov, result.flux_xy_Cov);
}

SdssShapeAlgorithm::ResultMapper SdssShapeAlgorithm::makeResultMapper(
    afw::table::Schema & schema,
    std::string const & name,
    Control const & ctrl
) {
    return ResultMapper(schema, name, FULL_COVARIANCE, SIGMA_ONLY, FULL_COVARIANCE, ctrl);
}

template <typename T>
SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(
    afw::image::MaskedImage<T> const & exposure,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & position,
    Control const & ctrl
) {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Not implemented"
    );
}

template <typename T>
SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(
    afw::image::Image<T> const & exposure,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & position,
    Control const & ctrl
) {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Not implemented"
    );
}

template <typename T>
SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure.getMaskedImage(), *inputs.footprint, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(      \
        afw::image::MaskedImage<T> const & exposure,                    \
        afw::detection::Footprint const & footprint,                    \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(      \
        afw::image::Image<T> const & exposure,                          \
        afw::detection::Footprint const & footprint,                    \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(               \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

/// @endcond
