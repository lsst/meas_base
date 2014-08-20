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
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/meas/base/PsfFlux.h"

namespace lsst { namespace meas { namespace base {

PsfFluxAlgorithm::ResultMapper PsfFluxAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl
) {
    return ResultMapper(schema, name, SIGMA_ONLY);
}

template <typename T>
 void PsfFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    afw::geom::Point2D const & position,
    Result & result,
    Control const & ctrl
) {
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[NO_PSF].doc,
            NO_PSF
        );
    }
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
    typedef afw::detection::Psf::Pixel PsfPixel;
    typedef typename afw::image::MaskedImage<T>::Variance::Pixel VarPixel;
    ndarray::EigenView<PsfPixel,1,1,Eigen::ArrayXpr> model(
        afw::detection::flattenArray(
            fitRegion,
            psfImage->getArray(),
            psfImage->getXY0()
        )
    );
    ndarray::EigenView<T,1,1,Eigen::ArrayXpr> data(
        afw::detection::flattenArray(
            fitRegion,
            exposure.getMaskedImage().getImage()->getArray(),
            exposure.getXY0()
        )
    );
    ndarray::EigenView<VarPixel,1,1,Eigen::ArrayXpr> variance(
        afw::detection::flattenArray(
            fitRegion,
            exposure.getMaskedImage().getVariance()->getArray(),
            exposure.getXY0()
        )
    );
    PsfPixel alpha = model.matrix().squaredNorm();
    result.flux = model.matrix().dot(data.matrix().template cast<PsfPixel>()) / alpha;
    // If we're not using per-pixel weights to compute the flux, we'll still want to compute the
    // variance as if we had, so we'll apply the weights to the model vector now, and update alpha.
    result.fluxSigma = std::sqrt(model.square().matrix().dot(variance.matrix().template cast<PsfPixel>()))
        / alpha;
    if (!utils::isfinite(result.flux) || !utils::isfinite(result.fluxSigma)) {
        throw LSST_EXCEPT(PixelValueError, "Invalid pixel value detected in image.");
    }
    if (result.getFlag(EDGE)) {
        throw LSST_EXCEPT(
            MeasurementError,
            getFlagDefinitions()[EDGE].doc,
            EDGE
        );
    }
}

template <typename T>
 void PsfFluxAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Result & result,
    Control const & ctrl
) {
    apply(exposure, inputs.position, result, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template  void PsfFluxAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        afw::geom::Point2D const & position,                            \
        Result & result,                                          \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
     void PsfFluxAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Result & result,                                          \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base

