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

#ifndef LSST_MEAS_BASE_SdssShape_h_INCLUDED
#define LSST_MEAS_BASE_SdssShape_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

class SdssShapeControl {
public:

    LSST_CONTROL_FIELD(background, double, "Additional value to add to background");
    LSST_CONTROL_FIELD(maxIter, int, "Maximum number of iterations");
    LSST_CONTROL_FIELD(maxShift, double, "Maximum centroid shift");
    LSST_CONTROL_FIELD(tol1, float, "Convergence tolerance for e1,e2");
    LSST_CONTROL_FIELD(tol2, float, "Convergence tolerance for FWHM");

    SdssShapeControl() : background(0.0), maxIter(100), maxShift(), tol1(1E-5), tol2(1E-4) {}

};

class SdssShapeAlgorithm {
public:

    enum FlagBits { N_FLAGS=0 };

    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions();

    typedef SdssShapeControl Control;
    typedef SimpleResult3<SdssShapeAlgorithm,ShapeResult,CentroidResult,FluxResult> Result;
    typedef SimpleResultMapper3<SdssShapeAlgorithm,ShapeResultMapper,CentroidResultMapper,FluxResultMapper>
                ResultMapper;
    typedef AlgorithmInput2 Input;

    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & name,
        Control const & ctrl=Control()
    );

    // Main implementation
    template <typename T>
    static Result apply(
        afw::image::MaskedImage<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    // Limited implementation - no variance plane, and hence no uncertainties
    template <typename T>
    static Result apply(
        afw::image::Image<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position,
        Control const & ctrl=Control()
    );

    // Forwarding overload called by plugin framework
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssShape_h_INCLUDED
