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
#include <iostream>
#include <cmath>
#include <numeric>

#include "ndarray/eigen.h"

#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/base/ApertureFlux.h"

namespace lsst { namespace meas { namespace base {

namespace {

template <typename MaskedImageT>
class FootprintFlux : public afw::detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : afw::detection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0.0), _sumVar(0.0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = _sumVar = 0.0;
    }
    void reset(afw::detection::Footprint const&) {}        

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

    // Return the Footprint's flux
    double getSum() const { return _sum; }

    // Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:
    double _sum;
    double _sumVar;
};

template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public afw::detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage,          // The image the source lives in
                        typename WeightImageT::Ptr wimage    // The weight image
                       ) : afw::detection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0.0), _sumVar(0.0), _x0(0), _y0(0) {}

    // Reset everything for a new Footprint
    void reset(afw::detection::Footprint const& foot) {
        _sum = _sumVar = 0.0;

        afw::geom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pex::exceptions::LengthError,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size for "
                                             "%d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }

    void reset() {}

    // method called for each pixel by apply()
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


} // anonymous

ApFluxComponentMapper::ApFluxComponentMapper(
    afw::table::Schema & schema,
    std::string const & prefix,
    UncertaintyEnum uncertainty
) {
    _flux = schema.addField(
        afw::table::Field<Flux>(prefix + "_flux", "measured flux", "dn"),
        true
    );
    if (uncertainty != NO_UNCERTAINTY) {
        _fluxSigma = schema.addField(
            afw::table::Field<ErrElement>(
                prefix + "_fluxSigma", "1-sigma uncertainty on " + prefix + "_flux", "dn"
            ),
            true
         );
    }
}

void ApFluxComponentMapper::apply(afw::table::BaseRecord & record, FluxComponent const & result) const {
    record.set(_flux, result.flux);
    if (_fluxSigma.isValid()) {
        record.set(_fluxSigma, result.fluxSigma);
    }
}

ApertureFluxExtras::ApertureFluxExtras(){};

ApertureFluxExtrasMapper::ApertureFluxExtrasMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ApertureFluxControl const & control
    ) :
    _nApertures(
        schema.addField(
            afw::table::Field<int>(
                prefix + "_nApertures", "number of apertures measured", ""
            ), true
        )
    )
{
    // Using control.radii, create a FluxComponentMapper for each radius
    _fluxComponentMapperVector.clear();
    for (unsigned int i = 0; i < control.radii.size(); i++)
    {
        std::stringstream sstm;
        sstm << prefix << '.' << i;
        FluxComponentMapper fluxComponentMapper(schema, sstm.str(), SIGMA_ONLY);
        _fluxComponentMapperVector.push_back(fluxComponentMapper);
    }
};

void ApertureFluxExtrasMapper::apply(afw::table::BaseRecord & record, ApertureFluxExtras const & result) const
{
    record.set(_nApertures, _fluxComponentMapperVector.size());
    for (unsigned int i = 0; i < _fluxComponentMapperVector.size(); i++)
    {
        _fluxComponentMapperVector[i].apply(record, *result.fluxComponentVector[i].get());
    }
};

/**
 * Create the object that controls aperture photometry
 */
ApertureFluxControl::ApertureFluxControl() : radii()
{
}

ApertureFluxAlgorithm::ResultMapper ApertureFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl)
{
    return ResultMapper(schema, name, ctrl);
}

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & position,
    Control const & ctrl
) {
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;

    // The following code plus the Footprint helper classes were taken from the meas_algorithms
    // ApertureFlux class with as few changes as possible, so that we can integrate later changes

    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;

    MaskedImageT const& mimage = exposure.getMaskedImage();

    double const xcen = position.getX();   ///< object's column position
    double const ycen = position.getY();   ///< object's row position

    int const ixcen = afw::image::positionToIndex(xcen);
    int const iycen = afw::image::positionToIndex(ycen);

    // BBox for data image
    lsst::afw::geom::BoxI imageBBox(mimage.getBBox(afw::image::PARENT));

    Result result;

    for (unsigned int i = 0; i <  ctrl.radii.size(); ++i) {
        FootprintFlux<MaskedImageT> fluxFunctor(mimage);
        afw::detection::Footprint const foot(afw::geom::PointI(ixcen, iycen), ctrl.radii[i], imageBBox);
        boost::shared_ptr<lsst::meas::base::FluxComponent> ptr
            = boost::shared_ptr<lsst::meas::base::FluxComponent>(new FluxComponent());
        result.fluxComponentVector.push_back(ptr);
        try {
            fluxFunctor.apply(foot);
            ptr.get()->flux = fluxFunctor.getSum();
            ptr.get()->fluxSigma = ::sqrt(fluxFunctor.getSumVar());
        }
        // this exception indicates that the footprint extends over the edge of the image
        // set the bit, but don't invalidate any successful measurements
        catch (lsst::pex::exceptions::LengthError &) {
            result.setFlag(EDGE);
        }
    }

    // end of meas_algorithms code

    return result;
}

template <typename T>
ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.position, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    ApertureFluxAlgorithm::Result ApertureFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

/// @endcond

