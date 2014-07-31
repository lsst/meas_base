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
#include <iostream>
#include <cmath>
#include <numeric>
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/algorithms/ApertureFluxTemplates.h"


// Doxygen gets confused and generates warnings when trying to map the definitions here to their
// declarations, but we want to put the source code itself in the HTML docs, so we just tell it
// not to look for any documentation comments here.
/// @cond SOURCE_FILE

namespace lsst { namespace meas { namespace base {


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

void ApFluxComponentMapper::apply(afw::table::BaseRecord & record, FluxComponent const & result) const { //changed pgee
    record.set(_flux, result.flux);
    if (_fluxSigma.isValid()) {
        record.set(_fluxSigma, result.fluxSigma);
    }
}

// pgee changed ApFluxComponent::ApFluxComponent() :
//    flux(std::numeric_limits<float>::quiet_NaN()),
//    fluxSigma(std::numeric_limits<float>::quiet_NaN())
//{}
ApertureFluxExtras::ApertureFluxExtras(){}; ///< Constructor; initializes everything to NaN

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
    // Using the could of radii on control.radii, create a FluxComponentMapper for each radius
    _fluxComponentMapperVector.clear();
    for (unsigned int i = 0; i < control.radii.size(); i++)
    {
        std::stringstream sstm;
        sstm << prefix << '.' << i;
        FluxComponentMapper fluxComponentMapper(schema, sstm.str(), SIGMA_ONLY);  // changed pgee
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
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_PSF].doc,
            NO_PSF
        );
    }
    Result result;
    PTR(afw::detection::Psf::Image) psfImage = psf->computeImage(position);
    afw::geom::Box2I fitBBox = psfImage->getBBox(afw::image::PARENT);
    fitBBox.clip(exposure.getBBox(afw::image::PARENT));
    if (fitBBox != psfImage->getBBox(afw::image::PARENT)) {
        result.setFlag(EDGE);
    }
    afw::detection::Footprint fitRegion(fitBBox);
    if (!ctrl.badMaskPlanes.empty()) {
        afw::image::MaskPixel badBits = 0x0;
        for (
            std::vector<std::string>::const_iterator i = ctrl.badMaskPlanes.begin();
            i != ctrl.badMaskPlanes.end();
            ++i
        ) {
            badBits |= exposure.getMaskedImage().getMask()->getPlaneBitMask(*i);
        }
        fitRegion.intersectMask(*exposure.getMaskedImage().getMask(), badBits);
    }
    if (fitRegion.getArea() == 0) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_GOOD_PIXELS].doc,
            NO_GOOD_PIXELS
        );
    }

    ndarray::EigenView<T,1,1,Eigen::ArrayXpr> data(
        afw::detection::flattenArray(
            fitRegion,
            exposure.getMaskedImage().getImage()->getArray(),
            exposure.getXY0()
        )
    );
   
    if (!utils::isfinite(data.sum())) {
        throw LSST_EXCEPT(PixelValueError, "Invalid pixel value detected in image.");
    }
    
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

    /* ******************************************************* */
    // Aperture flux
    result.fluxComponentVector.clear();

    for (unsigned int i = 0; i <  ctrl.radii.size(); ++i) {
        algorithms::FootprintFlux<MaskedImageT> fluxFunctor(mimage);
        afw::detection::Footprint const foot(afw::geom::PointI(ixcen, iycen), ctrl.radii[i], imageBBox);
        boost::shared_ptr<lsst::meas::base::FluxComponent> ptr = boost::shared_ptr<lsst::meas::base::FluxComponent>(new FluxComponent()); // changed pgee
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

