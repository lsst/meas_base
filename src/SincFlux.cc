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

#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/SincCoeffs.h"
#include "lsst/meas/base/SincFlux.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwGeom = lsst::afw::geom;

namespace lsst { namespace meas { namespace base {

namespace {

template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public afwDet::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(
        MaskedImageT const& mimage, ///< The image the source lives in
        PTR(WeightImageT) wimage    ///< The weight image
    ) :
        afwDet::FootprintFunctor<MaskedImageT>(mimage),
        _wimage(wimage),
        _sum(0.0), _sumVar(0.0),
        _x0(wimage->getX0()), _y0(wimage->getY0()) {}

    /// @brief Reset everything for a new Footprint
    void reset() {}
    void reset(afwDet::Footprint const& foot) {
        _sum = 0.0;
        _sumVar = 0.0;

        afwGeom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthError,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                             "for %d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
    ) {
        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        typename MaskedImageT::Image::Pixel vval = iloc.variance(0, 0);
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum    += wval*ival;
        _sumVar += wval*wval*vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    double getSumVar() const { return _sumVar; }

private:
    PTR(WeightImageT const) _wimage;                  // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;                                   // sum of the variance
    int _x0, _y0;                                     // the origin of the current Footprint
};

/************************************************************************************************************/

/**
 * Workhorse routine to calculate elliptical aperture fluxes
 */
template<typename MaskedImageT>
std::pair<double, double>
calculateSincApertureFlux(MaskedImageT const& mimage, afw::geom::ellipses::Ellipse const& ellipse,
                          double const innerFactor)
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;

    // BBox for data image
    afwGeom::BoxI imageBBox(mimage.getBBox());


    // make the coeff image
    // compute c_i as double integral over aperture def g_i(), and sinc()
    CONST_PTR(Image) cimage0 = SincCoeffs<Pixel>::get(ellipse.getCore(), innerFactor);

    // as long as we're asked for the same radius, we don't have to recompute cimage0
    // shift to center the aperture on the object being measured
    ImagePtr cimage = afwMath::offsetImage(*cimage0, ellipse.getCenter().getX(), ellipse.getCenter().getY());
    afwGeom::BoxI bbox(cimage->getBBox());
#if 0
    // I (Steve Bickerton) think this should work, but doesn't.
    // For the time being, I'll do the bounds check here
    // ... should determine why bbox/image behaviour not as expected.
    afwGeom::BoxI mbbox(mimage.getBBox());
    bbox.clip(mbbox);
    afwGeom::Point2I cimXy0(cimage->getXY0());
    bbox.shift(-cimage->getX0(), -cimage->getY0());
    cimage = typename Image::Ptr(new Image(*cimage, bbox, false));
    cimage->setXY0(cimXy0);
#else
    int x1 = (cimage->getX0() < mimage.getX0()) ? mimage.getX0() : cimage->getX0();
    int y1 = (cimage->getY0() < mimage.getY0()) ? mimage.getY0() : cimage->getY0();
    int x2 = (cimage->getX0() + cimage->getWidth() > mimage.getX0() + mimage.getWidth()) ?
        mimage.getX0() + mimage.getWidth() - 1 : cimage->getX0() + cimage->getWidth() - 1;
    int y2 = (cimage->getY0() + cimage->getHeight() > mimage.getY0() + mimage.getHeight()) ?
        mimage.getY0() + mimage.getHeight() - 1 : cimage->getY0() + cimage->getHeight() - 1;

    // if the dimensions changed, put the image in a smaller bbox
    if ( (x2 - x1 + 1 != cimage->getWidth()) || (y2 - y1 + 1 != cimage->getHeight()) ) {
        bbox = afwGeom::BoxI(afwGeom::Point2I(x1 - cimage->getX0(), y1 - cimage->getY0()),
                             afwGeom::Extent2I(x2 - x1 + 1, y2 - y1 + 1));
        cimage = ImagePtr(new Image(*cimage, bbox, afwImage::LOCAL, false));

        // shift back to correct place
        cimage = afwMath::offsetImage(*cimage, x1, y1);
		bbox = afwGeom::BoxI(afwGeom::Point2I(x1, y1),
							  afwGeom::Extent2I(x2-x1+1, y2-y1+1));
    }
#endif

    // pass the image and cimage into the wfluxFunctor to do the sum
    FootprintWeightFlux<MaskedImageT, Image> wfluxFunctor(mimage, cimage);

    afwDet::Footprint foot(bbox, imageBBox);
    wfluxFunctor.apply(foot);
    flux = wfluxFunctor.getSum();
    fluxErr = ::sqrt(wfluxFunctor.getSumVar());

    return std::make_pair(flux, fluxErr);
}

} // anonymous

SincFluxAlgorithm::SincFluxAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _fluxResultKey(
        FluxResultKey::addFields(schema, name, "flux from Sinc Flux algorithm")
    ),
    _centroidExtractor(schema, name)
{
    static boost::array<FlagDefinition,N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag, set if anything went wrong"},
    }};
    _flagHandler = FlagHandler::addFields(schema, name, flagDefs.begin(), flagDefs.end());
}

void SincFluxAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    // get the value from the centroid slot only
    afw::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    FluxResult result;

    afw::geom::ellipses::Axes const axes(_ctrl.radius2, _ctrl.radius2);
    std::pair<double, double> fluxes = calculateSincApertureFlux(exposure.getMaskedImage(),
                                                                 afw::geom::ellipses::Ellipse(axes, center),
                                                                 _ctrl.radius1/_ctrl.radius2);
    double flux = fluxes.first;
    double fluxErr = fluxes.second;
    result.flux = flux;
    result.fluxSigma = fluxErr;

    //  End of meas_algorithms code

    measRecord.set(_fluxResultKey, result);
    _flagHandler.setValue(measRecord, FAILURE, false);
}


void SincFluxAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}} // namespace lsst::meas::base

