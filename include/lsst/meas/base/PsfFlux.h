// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST.
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

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Transform.h"

namespace lsst {
namespace meas {
namespace base {

/**
 *  @brief A C++ control class to handle PsfFluxAlgorithm's configuration
 *
 *  In C++, we define Control objects to handle configuration information.  Using the LSST_CONTROL_FIELD
 *  macro and lsst.pex.config.wrap.makeConfigClass, we can turn these into more full-featured Config classes
 *  in Python.  While the user will usually interact with the Config class, the plugin wrapper system will
 *  turn Config instances into Control instances when passing them to C++.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class PsfFluxControl {
public:
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be excluded from the fit");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    PsfFluxControl() {}
};

/**
 *  @brief A measurement algorithm that estimates flux using a linear least-squares fit with the Psf model
 *
 *  The PsfFlux algorithm is extremely simple: we do a least-squares fit of the Psf model (evaluated
 *  at a given position) to the data.  For point sources, this provides the optimal flux measurement
 *  in the limit where the Psf model is correct.  We do not use per-pixel weights in the fit, as this
 *  results in bright stars being fit with a different effective profile than faint stairs.
 */
class PsfFluxAlgorithm : public SimpleAlgorithm {
public:
    // Structures and routines to manage flaghandler
    static FlagDefinitionList const& getFlagDefinitions();
    static FlagDefinition const FAILURE;
    static FlagDefinition const NO_GOOD_PIXELS;
    static FlagDefinition const EDGE;

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef PsfFluxControl Control;

    PsfFluxAlgorithm(Control const& ctrl, std::string const& name, afw::table::Schema& schema,
                     std::string const& logName = "");

    virtual void measure(afw::table::SourceRecord& measRecord,
                         afw::image::Exposure<float> const& exposure) const;

    virtual void fail(afw::table::SourceRecord& measRecord, MeasurementError* error = nullptr) const;

private:
    Control _ctrl;
    FluxResultKey _fluxResultKey;
    afw::table::Key<float> _areaKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
};

class PsfFluxTransform : public FluxTransform {
public:
    typedef PsfFluxControl Control;
    PsfFluxTransform(Control const& ctrl, std::string const& name, afw::table::SchemaMapper& mapper);
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_PsfFlux_h_INCLUDED
