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

#ifndef LSST_MEAS_BASE_GaussianFlux_h_INCLUDED
#define LSST_MEAS_BASE_GaussianFlux_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Transform.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle GaussianFluxAlgorithm's configuration
 */
class GaussianFluxControl {
public:
    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    GaussianFluxControl() :
        background(0.0)
    {}
};

/**
 *  @brief A measurement algorithm that estimates flux using an elliptical Gaussian weight.
 *
 *  This algorithm computes flux as the dot product of an elliptical Gaussian weight function
 *  with the image.  The size and ellipticity of the weight function are determined using the
 *  SdssShape algorithm, or retreived from a named field.
 */
class GaussianFluxAlgorithm : public SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    struct Flags;
    static FlagDefinition const & getDefinition(std::string name);
    static std::string const & getFlagName(std::size_t number);
    static std::size_t getFlagCount();

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef GaussianFluxControl Control;

    GaussianFluxAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

private:

    Control _ctrl;
    FluxResultKey _fluxResultKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
    SafeShapeExtractor _shapeExtractor;
};

class GaussianFluxTransform : public FluxTransform {
public:
    typedef GaussianFluxControl Control;
    GaussianFluxTransform(Control const & ctrl, std::string const & name, afw::table::SchemaMapper & mapper) :
                          FluxTransform{name, mapper} { }
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_GaussianFlux_h_INCLUDED
