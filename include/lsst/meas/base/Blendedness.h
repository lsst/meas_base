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

namespace lsst { namespace meas { namespace base {

class BlendednessControl {
public:
    LSST_CONTROL_FIELD(
        doOld, bool,
        "Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)"
    );

    LSST_CONTROL_FIELD(
        doFlux, bool,
        "Whether to compute quantities related to the Gaussian-weighted flux"
    );

    LSST_CONTROL_FIELD(
        doShape, bool,
        "Whether to compute quantities related to the Gaussian-weighted shape"
    );

    LSST_CONTROL_FIELD(
        nSigmaWeightMax, double,
        "Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)"
    );

    BlendednessControl() :
        doOld(true),
        doFlux(true),
        doShape(true),
        nSigmaWeightMax(3.0)
    {}

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
    struct Flags;
    static FlagDefinition const & getDefinition(std::string name);
    static std::string const & getFlagName(std::size_t number);
    static std::size_t getFlagCount();

    typedef BlendednessControl Control;

    BlendednessAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    void measureChildPixels(
        afw::image::MaskedImage<float> const & image,
        afw::table::SourceRecord & child
    ) const;

    void measureParentPixels(
        afw::image::MaskedImage<float> const & image,
        afw::table::SourceRecord & child
    ) const;

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const {}

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const {}


private:

    void _measureMoments(
        afw::image::MaskedImage<float> const & image,
        afw::table::SourceRecord & child,
        afw::table::Key<double> const & fluxRawKey,
        afw::table::Key<double> const & fluxAbsKey,
        ShapeResultKey const & _shapeRawKey,
        ShapeResultKey const & _shapeAbsKey
      ) const;

    Control const _ctrl;
    afw::table::Key<double> _old;
    afw::table::Key<double> _fluxRaw;
    afw::table::Key<double> _fluxChildRaw;
    afw::table::Key<double> _fluxParentRaw;
    afw::table::Key<double> _fluxAbs;
    afw::table::Key<double> _fluxChildAbs;
    afw::table::Key<double> _fluxParentAbs;
    ShapeResultKey _shapeChildRaw;
    ShapeResultKey _shapeParentRaw;
    ShapeResultKey _shapeChildAbs;
    ShapeResultKey _shapeParentAbs;
    FlagHandler _flagHandler;
};

}}} // namespace lsst::meas::base

#endif // LSST_MEAS_BASE_Blendedness_h_INCLUDED
