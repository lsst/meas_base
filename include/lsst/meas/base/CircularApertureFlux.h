// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#ifndef LSST_MEAS_BASE_CircularApertureFlux_h_INCLUDED
#define LSST_MEAS_BASE_CircularApertureFlux_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/ApertureFlux.h"

namespace lsst { namespace meas { namespace base {

class CircularApertureFluxAlgorithm : public SimpleAlgorithm , public ApertureFluxAlgorithm{
public:
    typedef ApertureFluxControl Control;
    CircularApertureFluxAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema, daf::base::PropertySet & metadata);

private:

    // These are private so they doesn't shadow the other overloads in base classes;
    // we can still call it via the public method on the base class.  We could have
    // used a using declaration instead, but Swig had trouble with that here.

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;
    SafeCentroidExtractor _centroidExtractor;
    FlagHandler _flagHandler;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_CircularApertureFlux_h_INCLUDED
