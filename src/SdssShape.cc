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

boost::array<FlagDef,SdssShapeAlgorithm::N_FLAGS> const & SdssShapeAlgorithm::getFlagDefinitions() {
    static boost::array<FlagDef,N_FLAGS> const flagDefs = {};
    return flagDefs;
}

SdssShapeAlgorithm::ResultMapper SdssShapeAlgorithm::makeResultMapper(
    afw::table::Schema & schema,
    std::string const & name,
    Control const & ctrl
) {
    return ResultMapper(schema, name, NO_UNCERTAINTY, NO_UNCERTAINTY, NO_UNCERTAINTY);
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
