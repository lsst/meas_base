// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2016 LSST/AURA
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

#ifndef LSST_MEAS_BASE_Blendedness_h_INCLUDED
#define LSST_MEAS_BASE_Blendedness_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Image.h"
#include "lsst/meas/base/Transform.h"
#include "lsst/meas/base/FlagHandler.h"

namespace lsst {
namespace meas {
namespace base {

class BlendednessControl {
public:
    LSST_CONTROL_FIELD(
            doOld, bool,
            "Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)");

    LSST_CONTROL_FIELD(doFlux, bool, "Whether to compute quantities related to the Gaussian-weighted flux");

    LSST_CONTROL_FIELD(doShape, bool, "Whether to compute quantities related to the Gaussian-weighted shape");

    LSST_CONTROL_FIELD(nSigmaWeightMax, double,
                       "Radius factor that sets the maximum extent of the weight function (and hence the "
                       "flux measurements)");

    BlendednessControl() : doOld(true), doFlux(true), doShape(true), nSigmaWeightMax(3.0) {}
};

/**
 *  Compute metrics that measure how blended objects are.
 *
 *  Blendedness is initialized once for a given Schema, then run repeatedly
 *  by calls to measureChildPixels and measureParentPixels (in any order, possibly
 *  with multiple sources interleaved).  Since it needs access to both the image with
 *  with noise and the noise replaced children it cannot use the standard plugin
 *  interface.
 */
class BlendednessAlgorithm : public SimpleAlgorithm {
public:
    // Structures and routines to manage flaghandler
    static FlagDefinitionList const& getFlagDefinitions();
    static FlagDefinition const FAILURE;
    static FlagDefinition const NO_CENTROID;
    static FlagDefinition const NO_SHAPE;

    typedef BlendednessControl Control;

    BlendednessAlgorithm(Control const& ctrl, std::string const& name, afw::table::Schema& schema);

    /**
     *  Compute the posterior expectation value of the true instrumental flux in a pixel
     *  from its (Gaussian) likelihood and a flat nonnegative prior.
     *
     *  This computes
     *  @f[
     *       \frac{\int_0^\infty f \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(f-z)^2}{2\sigma^2}} df}
                  {\int_0^\infty \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(f-z)^2}{2\sigma^2}} df}
     *  @f]
     *  where @f$z@f$ is the (noisy) pixel value and @f$\sigma^2@f$ is the pixel variance.  This
     *  approaches @f$z@f$ when @f$z \gg \sigma@f$ and @f$0@f$ when @f$z \ll -\sigma@f$.
     *
     *  We use single precision here for performance reasons; this function is called in a loop
     *  over single-precision pixels, and relies on a number of calls to exp and erfc, which are
     *  much faster in single precision.
     */
    static float computeAbsExpectation(float data, float variance);

    /**
     *  Compute the bias induced by using the absolute value of a pixel instead of its value.
     *
     *  The computation assumes the true distribution for the pixel is a Gaussian
     *  with mean mu and the given variance.  To compute mu from the data, use
     *  computeAbsExpectation.
     *
     *  This computes
     *  @f[
     *      \sqrt{\frac{2}{\pi}}\sigma e^{-\frac{\mu^2}{2\sigma^2}}
     *          - \mu\,\mathrm{erfc}\left(\frac{\mu}{\sqrt{2}\sigma}\right)
     *  @f]
     *  where @f$\mu@f$ is the mean of the underlying distribution and @f$\sigma^2@f$
     *  is its variance.
     *  See section 4.9.11 of Bosch, J. et al. 2018, PASJ, 70, S5 for further details.
     */
    static float computeAbsBias(float mu, float variance);

    void measureChildPixels(afw::image::MaskedImage<float> const& image,
                            afw::table::SourceRecord& child) const;

    void measureParentPixels(afw::image::MaskedImage<float> const& image,
                             afw::table::SourceRecord& child) const;

    virtual void measure(afw::table::SourceRecord& measRecord,
                         afw::image::Exposure<float> const& exposure) const {}

    virtual void fail(afw::table::SourceRecord& measRecord, MeasurementError* error = nullptr) const {}

private:
    void _measureMoments(afw::image::MaskedImage<float> const& image, afw::table::SourceRecord& child,
                         afw::table::Key<double> const& instFluxRawKey,
                         afw::table::Key<double> const& instFluxAbsKey, ShapeResultKey const& _shapeRawKey,
                         ShapeResultKey const& _shapeAbsKey) const;

    Control const _ctrl;
    afw::table::Key<double> _old;
    afw::table::Key<double> _raw;
    afw::table::Key<double> _instFluxChildRaw;
    afw::table::Key<double> _instFluxParentRaw;
    afw::table::Key<double> _abs;
    afw::table::Key<double> _instFluxChildAbs;
    afw::table::Key<double> _instFluxParentAbs;
    ShapeResultKey _shapeChildRaw;
    ShapeResultKey _shapeParentRaw;
    ShapeResultKey _shapeChildAbs;
    ShapeResultKey _shapeParentAbs;
    FlagHandler _flagHandler;
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // LSST_MEAS_BASE_Blendedness_h_INCLUDED
