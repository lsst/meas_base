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

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/base/NaiveFlux.h"

namespace lsst { namespace meas { namespace base {

namespace { //anonymous

template <typename MaskedImageT>
class FootprintFlux : public lsst::afw::detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : lsst::afw::detection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0.0), _sumVar(0.0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = _sumVar = 0.0;
    }
    void reset(lsst::afw::detection::Footprint const&) {}        

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int,                                   ///< column-position of pixel
                    int                                    ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = loc.image(0, 0);
        typename MaskedImageT::Variance::Pixel vval = loc.variance(0, 0);
        _sum += ival;
        _sumVar += vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

    /// Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:
    double _sum;
    double _sumVar;
};

template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public lsst::afw::detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage,          ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : lsst::afw::detection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0.0), _sumVar(0.0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset(lsst::afw::detection::Footprint const& foot) {
        _sum = _sumVar = 0.0;
        
        lsst::afw::geom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthError,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size for "
                                             "%d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    void reset() {}
    
    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        typename MaskedImageT::Variance::Pixel vval = iloc.variance(0, 0);
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
        _sumVar += wval*wval*vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the variance in the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;                                   // The variance of our desired sum
    int _x0, _y0;                                     // the origin of the current Footprint
};

            
/*****************************************************************************************************/
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
} // end anonymous namespace

NaiveFluxAlgorithm::ResultMapper NaiveFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
void NaiveFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & center,
    Result & result,
    Control const & ctrl
) {
    typename afw::image::Exposure<T>::MaskedImageT const& mimage = exposure.getMaskedImage();

    double const xcen = center.getX();   ///< object's column position
    double const ycen = center.getY();   ///< object's row position

    int const ixcen = afw::image::positionToIndex(xcen);
    int const iycen = afw::image::positionToIndex(ycen);

    // BBox for data image
    afw::geom::BoxI imageBBox(mimage.getBBox(afw::image::PARENT));

    /* ******************************************************* */
    // Aperture flux
    FootprintFlux<typename afw::image::Exposure<T>::MaskedImageT> fluxFunctor(mimage);
    afw::detection::Footprint const foot(
        afw::geom::PointI(ixcen, iycen),
        ctrl.radius,
        imageBBox
        );
    fluxFunctor.apply(foot);

    result.flux = fluxFunctor.getSum();
    result.fluxSigma = ::sqrt(fluxFunctor.getSumVar());
}

template <typename T>
void NaiveFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure, inputs.position, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void NaiveFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void NaiveFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base
