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
#include "lsst/meas/base/GaussianCentroid.h"
#include "lsst/meas/base/algorithms/all.h"


namespace lsst { namespace meas { namespace base {


GaussianCentroidAlgorithm::ResultMapper GaussianCentroidAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, NO_UNCERTAINTY);
}

template <typename T>
GaussianCentroidAlgorithm::Result GaussianCentroidAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Control const & ctrl
) {
    Result result;

    // This code has been moved essentially without change from meas_algorithms
    // The only changes were:
    // Change the exceptions to MeasurementErrors with the correct flag bits
    // Change to set values in result rather than in a source record.
    result.x = center.getX(); result.y = center.getY(); // better than NaN

    typedef afw::image::Image<T> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = static_cast<int>(center.getX() + 0.5);
    int y = static_cast<int>(center.getY() + 0.5);

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    FittedModel fit = twodg(image, x, y); // here's the fitter
    if (fit.params[FittedModel::PEAK] <= 0) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_PEAK].doc,
            NO_PEAK
        );
    }

    result.x = lsst::afw::image::indexToPosition(image.getX0()) + fit.params[FittedModel::X0];
    result.y = lsst::afw::image::indexToPosition(image.getY0()) + fit.params[FittedModel::Y0];

    // FIXME: should report uncertainty
    return result;
}

template <typename T>
GaussianCentroidAlgorithm::Result GaussianCentroidAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template GaussianCentroidAlgorithm::Result GaussianCentroidAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    GaussianCentroidAlgorithm::Result GaussianCentroidAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base


