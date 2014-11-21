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
#include "lsst/afw/table/Source.h"
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
PixelFlagsAlgorithm::PixelFlagsAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _centroidKey(
        CentroidResultKey::addFields(schema, name, "centroid from Naive Centroid algorithm", SIGMA_ONLY)
    ),
    _centroidExtractor(schema, name)
{
    static boost::array<FlagDefinition,N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag, set if anything went wrong"},
        {"flag_edge", "Could not use full PSF model image in fit because of proximity to exposure border"},
        {"flag_interpolated", "Interpolated pixel in the source footprint"},
        {"flag_interpolatedCenter", "Interpolated pixel in the source center"},
        {"flag_saturated", "Saturated pixel in the source footprint"},
        {"flag_saturatedCenter", "Saturated pixel in the source center"},
        {"flag_cr", "Cosmic ray in the source footprint"},
        {"flag_crCenter", "Cosmic ray in the source center"},
        {"flag_bad", "Bad pixel in the source footprint"}
    }};
    _flagHandler = FlagHandler::addFields(schema, name, flagDefs.begin(), flagDefs.end());
}  

void PixelFlagsAlgorithm::measure( 
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    
    afw::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    typedef typename afw::image::MaskedImage<float> MaskedImageT;
    MaskedImageT mimage = exposure.getMaskedImage();

    FootprintBits<MaskedImageT> func(mimage);

//  Catch NAN in centroid estimate
    if (lsst::utils::isnan(center.getX()) || lsst::utils::isnan(center.getY())) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Center point passed to PixelFlagsALgorithm is NaN");
    }
//  Catch centroids off the image
    if (!mimage.getBBox().contains(afw::geom::Point2I(center))) {
       _flagHandler.setValue(measRecord, EDGE, true);
    }
    // Check for bits set in the source's Footprint
    afw::detection::Footprint const & footprint(*measRecord.getFootprint());
    func.apply(footprint);
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
        _flagHandler.setValue(measRecord, EDGE, true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("BAD")) {
        _flagHandler.setValue(measRecord, BAD, true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        _flagHandler.setValue(measRecord, INTERPOLATED, true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        _flagHandler.setValue(measRecord, SATURATED, true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("CR")) {
        _flagHandler.setValue(measRecord, CR, true);
    }

    // Check for bits set in the 3x3 box around the center
    afw::geom::Point2I llc(afw::image::positionToIndex(center.getX())-1, afw::image::positionToIndex(center.getY()) - 1);

    afw::detection::Footprint const middle(afw::geom::BoxI(llc, afw::geom::ExtentI(3))); // central 3x3
    func.apply(middle);
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        _flagHandler.setValue(measRecord, INTERPOLATED_CENTER, true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        _flagHandler.setValue(measRecord, SATURATED_CENTER, true);
    }
    if (func.getBits() & MaskedImageT::Mask::getPlaneBitMask("CR")) {
        _flagHandler.setValue(measRecord, CR_CENTER, true);
    }
    _flagHandler.setValue(measRecord, FAILURE, false);
}

void PixelFlagsAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}} // namespace lsst::meas::base

