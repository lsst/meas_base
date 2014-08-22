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
#include "lsst/meas/base/NaiveCentroid.h"


namespace lsst { namespace meas { namespace base {


NaiveCentroidAlgorithm::ResultMapper NaiveCentroidAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, NO_UNCERTAINTY);
}

template <typename T>
void NaiveCentroidAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Result & result,
    Control const & ctrl
) {

    // This code has been moved essentially without change from meas_algorithms
    // The only changes were:
    // Change the exceptions to MeasurementErrors with the correct flag bits
    // Change to set values in result rather than in a source record.
    result.x = center.getX(); result.y = center.getY(); // better than NaN
    typedef afw::image::Image<T> ImageT;
    ImageT const& image = *exposure.getMaskedImage().getImage();

    int x = center.getX();  // FIXME: this is different from GaussianCentroid and SdssCentroid here,
    int y = center.getY();  //        and probably shouldn't be.

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {

        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[EDGE].doc,
            EDGE
        );
    }

    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1))
        - 9 * ctrl.background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_COUNTS].doc,
            NO_COUNTS
        );
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    result.x = lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x / sum;
    result.y = lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y / sum;

    // FIXME: should report uncertainty
}

template <typename T>
void NaiveCentroidAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure, inputs.position, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void NaiveCentroidAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void NaiveCentroidAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base


