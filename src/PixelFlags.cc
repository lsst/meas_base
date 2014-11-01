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
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/meas/base/PixelFlags.h"

namespace lsst { namespace meas { namespace base {
namespace {

template <typename MaskedImageT>
class FootprintBits : public afw::detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintBits(MaskedImageT const& mimage) :
        afw::detection::FootprintFunctor<MaskedImageT>(mimage), _bits(0)
    {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _bits = 0x0;
    }
    virtual void reset(afw::detection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        _bits |= loc.mask(0, 0);
    }

    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    typename MaskedImageT::Mask::Pixel _bits;
};
}  // end anonymous namespace

PixelFlagsAlgorithm::ResultMapper PixelFlagsAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
void PixelFlagsAlgorithm::apply(
    afw::image::MaskedImage<T> const & mimage,
    afw::geom::Point2D const & center,
    afw::detection::Footprint const & footprint,
    Result & result,
    Control const & ctrl
) {
    typedef typename afw::image::MaskedImage<T> MaskedImageT;
    FootprintBits<MaskedImageT> func(mimage);

//  Catch centroids off the image or NAN
    if (lsst::utils::isnan(center.getX()) || lsst::utils::isnan(center.getY())) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Center point passed to PixelFlagsALgorithm is NaN");
    } 
    if (!mimage.getBBox().contains(afw::geom::Point2I(center) -
                                                     afw::geom::Extent2I(mimage.getXY0()))) {
       result.setFlag(EDGE);
    }

    // Check for bits set in the source's Footprint
    func.apply(footprint);
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
        result.setFlag(EDGE);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("BAD")) {
        result.setFlag(BAD);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        result.setFlag(INTERPOLATED);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        result.setFlag(SATURATED);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("CR")) {
        result.setFlag(CR);
    }

    // Check for bits set in the 3x3 box around the center
    afw::geom::Point2I llc(afw::image::positionToIndex(center.getX())-1, afw::image::positionToIndex(center.getY()) - 1);

    afw::detection::Footprint const middle(afw::geom::BoxI(llc, afw::geom::ExtentI(3))); // central 3x3
    func.apply(middle);
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        result.setFlag(INTERPOLATED_CENTER);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        result.setFlag(SATURATED_CENTER);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("CR")) {
        result.setFlag(CR_CENTER);
    }
}

template <typename T>
void PixelFlagsAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure.getMaskedImage(), inputs.position, *inputs.footprint, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void PixelFlagsAlgorithm::apply(        \
        afw::image::MaskedImage<T> const & mimage,                       \
        afw::geom::Point2D const & position,                            \
        afw::detection::Footprint const & footprint,                    \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void PixelFlagsAlgorithm::apply(                 \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    );

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

