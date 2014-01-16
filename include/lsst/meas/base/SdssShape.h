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

// SdssShape actually measures a flux and a centroid as well, so our result struct includes all three
class SdssShapeAlgorithmResult : public ShapeAlgorithmResult,
                                 public CentroidAlgorithmResult,
                                 public FluxAlgorithmResult
{
public:
    SdssShapeAlgorithmResult() {}
};

class SdssShapeAlgorithmResultMapper : public ShapeAlgorithmMapper,
                                       public CentroidAlgorithmMapper,
                                       public FluxAlgorithmMapper
// n.b. Multiple inheritance here is ok even though there's a diamond and we didn't use virtual inheritance:
// BaseAlgorithmMapper is designed such that if you inherit it multiple times, and you hence get multiple
// "copies" of it, those copies all refer to the same Schema entry, so you just end up with a couple of
// extra lightweight Key data members that don't do anything: a small price to pay for the boilerplate
// reduction we get by using multiple inheritance here.
// Also, sometimes we might want to get multiple copies of the base class (if we pass different prefixes
// to different constructors); that works fine too.
{
public:

    explicit SdssShapeAlgorithmResultMapper(afw::table::Schema & schema, std::string const & name);

    // Transfer values from the result struct to the record, and clear the failure flag field.
    void apply(afw::table::BaseRecord & record, SdssShapeAlgorithmResult const & result);

    // Set the failure flag field.
    void fail(afw::table::BaseRecord & record);

private:
    // TODO: more flags
};

class SdssShapeAlgorithm {
public:

    typedef SdssShapeControl Control;
    typedef SdssShapeAlgorithmResult Result;
    typedef SdssShapeAlgorithmResultMapper ResultMapper;
    typedef AlgorithmInput2 Input;

    static bool const HAS_CONTROL = true;

    static ResultMapper makeResultMapper(
        Control const & ctrl,
        afw::table::Schema & schema,
        std::string const & name
    );

    // Main implementation
    template <typename T>
    static Result apply(
        Control const & ctrl,
        afw::image::MaskedImage<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position
    );

    // Limited implementation - no variance plane, and hence no uncertainties
    template <typename T>
    static Result apply(
        Control const & ctrl,
        afw::image::Image<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position
    );

    // Forwarding overload called by plugin framework
    template <typename T>
    static Result apply(
        Control const & ctrl,
        afw::image::Exposure<T> const & exposure,
        Input const & inputs
    ) {
        return apply(ctrl, exposure.getMaskedImage(), *inputs.footprint, inputs.position);
    }

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_SdssShape_h_INCLUDED
