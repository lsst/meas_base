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
#include <numeric>
#include <cmath>
#include <functional>

#include "lsst/meas/base/PsfFlux.h"
#include "lsst/meas/base/Results.h"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintFunctor.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace lsst { namespace meas { namespace base {

PsfFluxAlgorithm::ResultMapper PsfFluxAlgorithm::makeResultMapper(afw::table::Schema & schema) {
    return ResultMapper(schema, "flux.psf");
}
template <typename TargetImageT, typename WeightImageT>

class FootprintWeightFlux : public afwDetection::FootprintFunctor<TargetImageT> {
public:
    FootprintWeightFlux(TargetImageT const& mimage, ///< The image the source lives in
                        PTR(WeightImageT) wimage    ///< The weight image
                       ) : afwDetection::FootprintFunctor<TargetImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _sumVar(0), _x0(0), _y0(0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {}
    void reset(afwDetection::Footprint const& foot) {
        _sumVar = _sum = 0.0;

        afwGeom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                             "for %d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
   
    /// @brief method called for each pixel by apply()
    virtual void operator()(typename TargetImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        _callImpl(iloc, x, y, typename TargetImageT::image_category());
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

    /// Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:

    template <typename LocatorT>
    void _callImpl(LocatorT iloc, int x, int y, afw::image::detail::MaskedImage_tag) {
        double ival = iloc.image(0, 0);
        double vval = iloc.variance(0, 0);
        double wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
        _sumVar += wval*wval*vval;
    }

    template <typename LocatorT>
    void _callImpl(LocatorT iloc, int x, int y, afw::image::detail::Image_tag) {
        double ival = *iloc;
        double wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval * ival;
    }

    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;
    int _x0, _y0;                                     // the origin of the current Footprint
};
/**
 * Accumulate sum(x) and sum(x**2)
 */
template<typename T>
struct getSum2 {
    getSum2() : sum(0.0), sum2(0.0) {}
   
    getSum2& operator+(T x) {
        sum += x;
        sum2 += x*x;

        return *this;
    }
   
    double sum;                         // \sum_i(x_i)
    double sum2;                        // \sum_i(x_i^2)
};

template <typename TargetImageT>
std::pair<double,double> computePsfFlux(
    TargetImageT const image,
    PTR(afw::detection::Psf::Image) const & wimage,
    afw::geom::Point2D const & center
) {
    afwGeom::BoxI imageBBox(image.getBBox(afwImage::PARENT));
    FootprintWeightFlux<TargetImageT, afwDetection::Psf::Image> wfluxFunctor(image, wimage);
    // Build a rectangular Footprint corresponding to wimage
    afwDetection::Footprint foot(wimage->getBBox(afwImage::PARENT), imageBBox);
    wfluxFunctor.apply(foot);

    getSum2<afwDetection::Psf::Pixel> sum;
    sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);

    double flux = wfluxFunctor.getSum()*sum.sum/sum.sum2;
    double fluxErr = ::sqrt(wfluxFunctor.getSumVar())*::fabs(sum.sum)/sum.sum2;
    return std::make_pair(flux, fluxErr);
}


template <typename T>
PsfFluxAlgorithm::Result PsfFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::detection::Footprint const & footprint,
    afw::geom::Point2D const & center
) {

    // the following code was extracted from PsfFlux.cc in meas_algorithms
    CONST_PTR(afwDetection::Psf) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "No PSF provided for PSF photometry");
    }

    PTR(afwDetection::Psf::Image) psfImage;
    try {
        psfImage = psf->computeImage(center);
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)")
                            % center.getX() % center.getY()).str());
        throw e;
    }

    std::pair<double,double> pair = computePsfFlux(exposure.getMaskedImage(), psfImage, center);
    // end extracted code
    Result result = Result(pair.first, pair.second);
    return result;
}

template <typename T>
std::vector<PsfFluxAlgorithm::Result> PsfFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    std::vector<Input> const & inputs
) {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Not implemented"
    );
}

#define INSTANTIATE(T)                                                  \
    template PsfFluxAlgorithm::Result PsfFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::detection::Footprint const & footprint,                    \
        afw::geom::Point2D const & position                             \
    );                                                                  \
    template                                                            \
    PsfFluxAlgorithm::Result PsfFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs                                            \
    );                                                                  \
    template                                                            \
    std::vector<PsfFluxAlgorithm::Result> PsfFluxAlgorithm::apply(      \
        afw::image::Exposure<T> const & exposure,                       \
        std::vector<Input> const & inputs                               \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base
