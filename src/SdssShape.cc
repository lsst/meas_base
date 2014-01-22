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
#include "lsst/meas/algorithms/detail/SdssShape.h"

namespace detail=lsst::meas::algorithms::detail;

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
    Control const & control,
    afw::image::MaskedImage<T> const & mimage,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & center
) {

    Result result = SdssShapeAlgorithmResult();
    //xsource.set(getKeys().flag, true); // say we've failed so that's the result if we throw
    //xsource.set(_centroidKeys.flag, true); // say we've failed so that's the result if we throw
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;

    double xcen = center.getX();         // object's column position
    double ycen = center.getY();         // object's row position

    xcen -= mimage.getX0();             // work in image Pixel coordinates
    ycen -= mimage.getY0();
    
    float shiftmax = 1;                 // Max allowed centroid shift \todo XXX set shiftmax from Policy
    if (shiftmax < 2) {
        shiftmax = 2;
    } else if (shiftmax > 10) {
        shiftmax = 10;
    }

    detail::SdssShapeImpl shapeImpl;
    bool anyFlags = false;
    try {
        (void)detail::getAdaptiveMoments(mimage, control.background, xcen, ycen, shiftmax, &shapeImpl,
                                         control.maxIter, control.tol1, control.tol2);
    } catch (pex::exceptions::Exception & err) {
        for (int n = 0; n < detail::SdssShapeImpl::N_FLAGS; ++n) {
             //source.set(_flagKeys[n], shapeImpl.getFlag(detail::SdssShapeImpl::Flag(n)));
             //anyFlags = anyFlags || source.get(_flagKeys[n]);
        }
        throw;
    }
/*
 * We need to measure the PSF's moments even if we failed on the object
 * N.b. This isn't yet implemented (but the code's available from SDSS)
 */
    for (int n = 0; n < detail::SdssShapeImpl::N_FLAGS; ++n) {
        //x source.set(_flagKeys[n], shapeImpl.getFlag(detail::SdssShapeImpl::Flag(n)));
        //x anyFlags = anyFlags || source.get(_flagKeys[n]);
    }
    
    result.centroid.value = lsst::afw::geom::Point2D(shapeImpl.getX(), shapeImpl.getY());
    result.centroid.cov(0,0) = shapeImpl.getXErr();
    result.centroid.cov(1,1) = shapeImpl.getYErr();
    result.centroid.cov(0,1) = 0,0;
    result.centroid.cov(1,0) = 0,0;
    // Jim:  So I need to save the fluxScale?  or multiply by?
    result.flux.value = shapeImpl.getI0();
    result.flux.err = shapeImpl.getI0Err();
    //x source.set(_centroidKeys.flag, anyFlags);
    result.value = afw::geom::ellipses::Quadrupole(shapeImpl.getIxx(), shapeImpl.getIyy(), shapeImpl.getIxy());
    // FIXME: should do off-diagonal covariance elements too
    result.cov(0,0) = shapeImpl.getIxxErr() * shapeImpl.getIxxErr();
    result.cov(1,1) = shapeImpl.getIyyErr() * shapeImpl.getIyyErr();
    result.cov(2,2) = shapeImpl.getIxyErr() * shapeImpl.getIxyErr();
    //xsource.set(getKeys().flag, anyFlags);
    return result;

}

template <typename T>
SdssShapeAlgorithm::Result SdssShapeAlgorithm::apply(
    Control const & ctrl,
    afw::image::Image<T> const & image,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & center
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
