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

#ifndef LSST_MEAS_BASE_PsfFlux_h_INCLUDED
#define LSST_MEAS_BASE_PsfFlux_h_INCLUDED

#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

class PsfFluxAlgorithm {
public:

    enum FlagBits { N_FLAGS=0 };

    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions();

    typedef SimpleResult1<PsfFluxAlgorithm,FluxResult> Result;
    typedef SimpleResultMapper1<PsfFluxAlgorithm,FluxResultMapper> ResultMapper;
    typedef AlgorithmInput2 Input; // type passed to apply in addition to Exposure.
    typedef NullControl Control; // control object - we don't have any configuration, so we use NullControl

    // Create an object that transfers Result values to a record associated with the given schema.
    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        Control const & ctrl=Control()
    );

    // On the apply() methods: I think it makes sense to template them, and do explicit instantation
    // for both double and float.  But the plugin framework will only use the float versions of apply.
    // Russell has suggested having different names for the different apply() methods, as overloading
    // doesn't work so well with Python, and I think that's a good idea.  But I don't know what to
    // call them all yet so I've just called them apply() (a new naming proposal will be sent to the
    // list shortly).

    // Method for single-object, single-frame measurement; this is where the real implementation goes
    // for the single-object, single-frame case.  This is the main one we expect non-framework users
    // to call.  Because PsfFlux needs to get the Psf from the Exposure to do its job, we don't have
    // a version that just takes a MaskedImage or Image.  For other algorithms we might, but the ones
    // the plugin framework call will always take an Exposure.
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        afw::detection::Footprint const & footprint,
        afw::geom::Point2D const & position
    );

    // This is the version that will be called by both the SFM framework and the forced measurement
    // framework, in single-object mode.  It will delegate to the above.
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );

    // This is the version that will be called by both the SFM framework and the forced measurement
    // framework, in multi-object mode.  It will have its own implementation, as that will be distinct
    // from the single-object case (well, we could actually implement the single-object case using
    // the multi-object code, but it's easy enough to write the single-object version in this case
    // that I won't).
    template <typename T>
    static std::vector<Result> applyN(
        afw::image::Exposure<T> const & exposure,
        std::vector<Input> const & inputs,
        Control const & ctrl=Control()
    );

    // I'm not gonna add apply() methods for multifit just yet.  They're trickier in multifit-specific
    // ways that I don't think will affect the overall plugin-wrapper-layer design.

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_PsfFlux_h_INCLUDED
