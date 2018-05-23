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

#include <cctype> // ::tolower
#include <algorithm> // std::transform
#include <cmath>

#include "ndarray/eigen.h"

#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/meas/base/PixelFlags.h"

namespace lsst { namespace meas { namespace base {
namespace {
template <typename MaskedImageT>
class FootprintBits {
public:
    explicit FootprintBits() : _bits(0) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _bits = 0x0;
    }

    void operator()(geom::Point2I const & point,
                    typename MaskedImageT::Mask::Pixel const & value) {
        _bits |= value;
    }

    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    typename MaskedImageT::Mask::Pixel _bits;
};

typedef afw::image::MaskedImage<float> MaskedImageF;

void updateFlags(PixelFlagsAlgorithm::KeyMap const & maskFlagToPixelFlag, const FootprintBits<MaskedImageF> & func,
                 afw::table::SourceRecord & measRecord) {
        for (auto const & i: maskFlagToPixelFlag) {
            try {
                if (func.getBits() & MaskedImageF::Mask::getPlaneBitMask(i.first)) {
                    measRecord.set(i.second, true);
                }
            }
            catch (pex::exceptions::InvalidParameterError & err) {
                throw LSST_EXCEPT(FatalAlgorithmError, err.what());
            }
        }
    }
} // end anonymous namespace

PixelFlagsAlgorithm::PixelFlagsAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
) : _ctrl(ctrl) {
    // Add generic keys first, which don't correspond to specific mask planes
    _generalFailureKey = schema.addField<afw::table::Flag>(name + "_flag",
                                        "general failure flag, set if anything went wring");
    _offImageKey = schema.addField<afw::table::Flag>(name+"_flag" + "_offimage",
                                        "Source center is off image");
    // Set all the flags that correspond to mask planes anywhere in the footprint
    _anyKeys["EDGE"] = schema.addField<afw::table::Flag>(name + "_flag_edge",
                                        "Source is outside usable exposure region (masked EDGE or NO_DATA)");
    _anyKeys["INTRP"] = schema.addField<afw::table::Flag>(name + "_flag_interpolated",
                                        "Interpolated pixel in the Source footprint");
    _anyKeys["SAT"] = schema.addField<afw::table::Flag>(name + "_flag_saturated",
                                        "Saturated pixel in the Source footprint");
    _anyKeys["CR"] = schema.addField<afw::table::Flag>(name + "_flag_cr",
                                        "Cosmic ray in the Source footprint");
    _anyKeys["BAD"] = schema.addField<afw::table::Flag>(name + "_flag_bad",
                                        "Bad pixel in the Source footprint");
    _anyKeys["SUSPECT"] = schema.addField<afw::table::Flag>(name + "_flag_suspect",
                                        "Source's footprint includes suspect pixels");
    // Flags that correspond to mask bits which occur in the center of the object
    _centerKeys["INTRP"] = schema.addField<afw::table::Flag>(name + "_flag_interpolatedCenter",
                                        "Interpolated pixel in the Source center");
    _centerKeys["SAT"] = schema.addField<afw::table::Flag>(name + "_flag_saturatedCenter",
                                        "Saturated pixel in the Source center");
    _centerKeys["CR"] = schema.addField<afw::table::Flag>(name + "_flag_crCenter",
                                        "Cosmic ray in the Source center");
    _centerKeys["SUSPECT"] = schema.addField<afw::table::Flag>(name + "_flag_suspectCenter",
                                        "Source's center is close to suspect pixels");

    // Read in the flags passed from the configuration, and add them to the schema
    for (auto const & i: _ctrl.masksFpCenter) {
        std::string maskName(i);
        std::transform(maskName.begin(), maskName.end(), maskName.begin(), ::tolower);
        _centerKeys[i] = schema.addField<afw::table::Flag>(name + "_flag_" + maskName + "Center",
                                        "Source center is close to "+ i + " pixels");
    }

    for (auto const & i: _ctrl.masksFpAnywhere) {
        std::string maskName(i);
        std::transform(maskName.begin(), maskName.end(), maskName.begin(), ::tolower);
        _anyKeys[i] = schema.addField<afw::table::Flag>(name + "_flag_" + maskName,
                                        "Source footprint includes " + i + " pixels");
    }
}

void PixelFlagsAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {

    MaskedImageF mimage = exposure.getMaskedImage();
    FootprintBits<MaskedImageF> func;

    // Check if the measRecord has a valid centroid key, i.e. it was centroided
    geom::Point2D center;
    if (measRecord.getTable()->getCentroidKey().isValid()) {
        center = measRecord.getCentroid();
        //  Catch NAN in centroid estimate
        if (std::isnan(center.getX()) || std::isnan(center.getY())) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                              "Center point passed to PixelFlagsAlgorithm is NaN");
        }
    }
    else {
        // Set the general failure flag because using the Peak might affect
        // the current measurement
        measRecord.set(_generalFailureKey, true);
        // Attempting to set pixel flags with no centroider, the center will
        // be determined though the peak pixel. The first peak in the
        // footprint (supplied by the measurement framework) should be the
        // highest peak so that one will be used as a proxy for the central
        // tendency of the distribution of flux for the record.
        PTR(afw::detection::Footprint) footprint = measRecord.getFootprint();
        // If there is no footprint or the footprint contains no peaks, throw
        // a runtime error.
        if (!footprint || footprint->getPeaks().empty()) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeError, "No footprint, or no footprint peaks detected");
        }
        else {
            center.setX(footprint->getPeaks().front().getFx());
            center.setY(footprint->getPeaks().front().getFy());
        }
    }

    //  Catch centroids off the image
    if (!mimage.getBBox().contains(geom::Point2I(center))) {
        measRecord.set(_offImageKey, true);
        measRecord.set(_anyKeys.at("EDGE"), true);
    }

    // Check for bits set in the source's Footprint
    afw::detection::Footprint const & footprint(*measRecord.getFootprint());
    footprint.getSpans()->clippedTo(mimage.getBBox())->applyFunctor(func, *(mimage.getMask()));

    // Set the EDGE flag if the bitmask has NO_DATA set
    try {
        if (func.getBits() & MaskedImageF::Mask::getPlaneBitMask("NO_DATA")) {
            measRecord.set(_anyKeys.at("EDGE"), true);
        }
    } catch (pex::exceptions::InvalidParameterError & err) {
        throw LSST_EXCEPT(FatalAlgorithmError, err.what());
    }

    // update the source record for the any keys
    updateFlags(_anyKeys, func, measRecord);

    // Check for bits set in the 3x3 box around the center
    geom::Point2I llc(afw::image::positionToIndex(center.getX()) - 1,
                           afw::image::positionToIndex(center.getY()) - 1);

    func.reset();
    auto spans = std::make_shared<afw::geom::SpanSet>(geom::Box2I(llc, geom::ExtentI(3)));
    afw::detection::Footprint const middle(spans); // central 3x3
    middle.getSpans()->clippedTo(mimage.getBBox())->applyFunctor(func, *(mimage.getMask()));

    // Update the flags which have to do with the center of the footprint
    updateFlags(_centerKeys, func, measRecord);

}

void PixelFlagsAlgorithm::fail(afw::table::SourceRecord & measRecord, MeasurementError * error) const {
    measRecord.set(_generalFailureKey, true);
}

}}} // namespace lsst::meas::base
