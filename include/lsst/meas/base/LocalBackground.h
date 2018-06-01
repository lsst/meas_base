// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018 AURA/LSST.
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

#ifndef LSST_MEAS_BASE_LocalBackground_h_INCLUDED
#define LSST_MEAS_BASE_LocalBackground_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Transform.h"

namespace lsst {
namespace meas {
namespace base {

/// Configuration of LocalBackgroundAlgorithm
class LocalBackgroundControl {
public:
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be excluded from the measurement");
    LSST_CONTROL_FIELD(annulusInner, float,
                       "Inner radius for background annulus as a multiple of the PSF sigma");
    LSST_CONTROL_FIELD(annulusOuter, float,
                       "Outer radius for background annulus as a multiple of the PSF sigma");
    LSST_CONTROL_FIELD(bgRej, float,
                       "Rejection threshold (in standard deviations) for background measurement");
    LSST_CONTROL_FIELD(bgIter, int, "Number of sigma-clipping iterations for background measurement");

    LocalBackgroundControl()
            :  // deliberately omitting DETECTED: objects shouldn't be present if we're deblending, and by
               // including DETECTED we wouldn't get any pixels for the background if we're on something
               // extended and over threshold which is exactly the case when we do want to measure the
               // background.
              badMaskPlanes({"BAD", "SAT", "NO_DATA"}),
              annulusInner(7.0),
              annulusOuter(15.0),
              bgRej(3.0),
              bgIter(3) {}
};

/// A measurement algorithm that estimates the local background value per pixel
class LocalBackgroundAlgorithm : public SimpleAlgorithm {
public:
    static FlagDefinitionList const& getFlagDefinitions();
    static FlagDefinition const FAILURE;
    static FlagDefinition const NO_GOOD_PIXELS;
    static FlagDefinition const NO_PSF;

    typedef LocalBackgroundControl Control;

    LocalBackgroundAlgorithm(Control const& ctrl, std::string const& name, afw::table::Schema& schema,
                             std::string const& logName = "");

    virtual void measure(afw::table::SourceRecord& measRecord,
                         afw::image::Exposure<float> const& exposure) const;

    virtual void fail(afw::table::SourceRecord& measRecord, MeasurementError* error = nullptr) const;

private:
    Control _ctrl;
    FluxResultKey _resultKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
    afw::math::StatisticsControl _stats;
};

class LocalBackgroundTransform : public FluxTransform {
public:
    typedef LocalBackgroundControl Control;
    LocalBackgroundTransform(Control const& ctrl, std::string const& name, afw::table::SchemaMapper& mapper);
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_LocalBackground_h_INCLUDED
