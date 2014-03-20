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
#include "ndarray/eigen.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/meas/base/SkyCoord.h"


// Doxygen gets confused and generates warnings when trying to map the definitions here to their
// declarations, but we want to put the source code itself in the HTML docs, so we just tell it
// not to look for any documentation comments here.
/// @cond SOURCE_FILE

namespace lsst { namespace meas { namespace base {


SkyCoordExtras::SkyCoordExtras() {}; ///< Constructor; initializes everything to NaN

SkyCoordExtrasMapper::SkyCoordExtrasMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        SkyCoordControl const & control
    )
// No field needed as we will be using the default coord
{};

// This has to be calculated during apply, because it is the only time we have the source record
void SkyCoordExtrasMapper::apply(afw::table::BaseRecord & record, SkyCoordExtras const & result) const
{
    afw::coord::Coord coord = record.get(_coord);
    //coord = *wcs.pixelToSky(getCentroid()));
    record.set(_coord, coord);
};

SkyCoordAlgorithm::ResultMapper SkyCoordAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl)
{
    return ResultMapper(schema, name, ctrl);
}

template <typename T>
SkyCoordAlgorithm::Result SkyCoordAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & position,
    Control const & ctrl
) {
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    PTR(afw::image::Wcs const) wcs = exposure.getWcs();
    if (!wcs) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_WCS].doc,
            NO_WCS
        );
    }
    Result result;

    // end of meas_algorithms code

    return result;
}

template <typename T>
SkyCoordAlgorithm::Result SkyCoordAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template SkyCoordAlgorithm::Result SkyCoordAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    SkyCoordAlgorithm::Result SkyCoordAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

/// @endcond

