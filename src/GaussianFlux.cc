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

#include "Eigen/Core"
#include "Eigen/LU"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/meas/base/GaussianFlux.h"
#include "lsst/meas/base/detail/SdssShapeImpl.h"

namespace lsst { namespace meas { namespace base {


/************************************************************************************************************/

GaussianFluxAlgorithm::ResultMapper GaussianFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
void GaussianFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & centroid,
    afw::geom::ellipses::Quadrupole const & shape,
    Result & result,
    Control const & ctrl
) {
    //  This code came straight out of the GaussianFlux.apply() in meas_algorithms with few changes
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    typename afw::image::Exposure<T>::MaskedImageT const& mimage = exposure.getMaskedImage();

    double const xcen = centroid.getX() - mimage.getX0(); ///< column position in image pixel coords
    double const ycen = centroid.getY() - mimage.getY0(); ///< row position


    detail::SdssShapeImpl sdss(centroid, shape);
    std::pair<double, double> fluxResult
        = detail::getFixedMomentsFlux(mimage, ctrl.background, xcen, ycen, sdss);
    result.flux =  fluxResult.first;
    result.fluxSigma = fluxResult.second;

    //  End of meas_algorithms code
}

template <typename T>
void GaussianFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure, inputs.position, inputs.shape, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void GaussianFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & centroid, \
        afw::geom::ellipses::Quadrupole const & shape, \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void GaussianFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

