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

#include "lsst/meas/base/SdssShape.h"

namespace lsst { namespace meas { namespace base {

SdssShapeAlgorithmResultMapper::SdssShapeAlgorithmResultMapper(afw::table::Schema & schema) :
    ShapeAlgorithmMapper(schema, "shape.sdss"),
    _centroid(schema, "shape.sdss.centroid"),
    _flux(schema, "shape.sdss.flux")
{}

void SdssShapeAlgorithmResultMapper::apply(
    afw::table::BaseRecord & record,
    SdssShapeAlgorithmResult const & result
) {
    ShapeAlgorithmMapper::apply(record, result);
    _centroid.apply(record, result.centroid);
    _flux.apply(record, result.flux);
}

void SdssShapeAlgorithmResultMapper::fail(afw::table::BaseRecord & record) {
    ShapeAlgorithmMapper::fail(record);
    _centroid.fail(record);
    _flux.fail(record);
}

SdssShapeAlgorithm::ResultMapper SdssShapeAlgorithm::makeResultMapper(afw::table::Schema & schema) {
    return ResultMapper(schema);
}

template <typename T>
SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(
    Control const & ctrl,
    afw::image::MaskedImage<T> const & exposure,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & position
) {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Not implemented"
    );
}

template <typename T>
SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(
    Control const & ctrl,
    afw::image::Image<T> const & exposure,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & position
) {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Not implemented"
    );
}

#define INSTANTIATE(T)                                                  \
    template SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(      \
        Control const & ctrl,                                           \
        afw::image::MaskedImage<T> const & exposure,                    \
        afw::detection::Footprint const & footprint,                    \
        afw::geom::Point2D const & position                             \
    );                                                                  \
    template SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(      \
        Control const & ctrl,                                           \
        afw::image::Image<T> const & exposure,                          \
        afw::detection::Footprint const & footprint,                    \
        afw::geom::Point2D const & position                             \
    );                                                                  \
    template                                                            \
    SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(               \
        Control const & ctrl,                                           \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base
