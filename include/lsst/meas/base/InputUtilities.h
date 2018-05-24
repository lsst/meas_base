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

#ifndef LSST_MEAS_BASE_InputUtilities_h_INCLUDED
#define LSST_MEAS_BASE_InputUtilities_h_INCLUDED

#include "lsst/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/meas/base/FlagHandler.h"

namespace lsst {
namespace meas {
namespace base {

/**
 *  Utility class for measurement algorithms that extracts a position from the Centroid slot and handles
 *  errors in a safe and consistent way.
 */
class SafeCentroidExtractor {
public:
    /**
     *  Construct the extractor, creating a flag alias that indicates failure in the input centroid
     *  by linking to the slot centroid flag.
     *
     *  @param[out] schema   Schema to which the alias should be added.  The "slot_Centroid" alias
     *                       must already be present in the Schema's AliasMap.
     *  @param[in]  name     The name of the algorithm; the flag alias added will be
     *                       "<name>_flag_badCentroid", or "<name>_flag_badInitialCentroid"
     *                       if isCentroider=true.
     *  @param[in]  isCentroider    Indicates whether the calling algorithm is itself a centroid
     *                              measurement algorithm.  If true,, falling back to the Peak
     *                              because there was no previous centroider or a previous centroider
     *                              failed will not cause the general failure flag of the current
     *                              algorithm to be set.
     */
    SafeCentroidExtractor(afw::table::Schema& schema, std::string const& name, bool isCentroider = false);

    /**
     *  Extract a position from the given record.
     *
     *  We use the Centroid slot if it is not NaN, and fall back to the Peak on the Footprint otherwise.
     *
     *  If the Centroid slot is not defined, we throw FatalAlgorithmError, as this indicates a configuration
     *  problem.
     *
     *  If the Centroid slot value is NaN and is *not* flagged (or there is no Centroid slot flag), we throw
     *  RuntimeError, which should cause the measurement framework to log a warning and set the current
     *  algorithm's general failure flag if it is allowed to propagate out of the algorithm implementation.
     *
     *  If the Centroid slot is NaN and there is no Footprint or Peak (even if the centroid is flagged),
     *  we also throw RuntimeError, as this indicates something wrong in another stage of the pipeline
     *  that should be addressed before measurement is run.
     *
     *  If the Centroid slot is flagged and we nevertheless obtain a usable position (either from a
     *  the Centroid slot itself or from a successful fall-back to the Peak), we set the current
     *  algorithm's general failure flag, but return the position as well, allowing it to continue while
     *  indicating that the result may not be reliable.
     */
    geom::Point2D operator()(afw::table::SourceRecord& record, FlagHandler const& flags) const;

private:
    std::string _name;
    bool _isCentroider;
};

/**
 *  Utility class for measurement algorithms that extracts an ellipse from the Shape slot and handles
 *  errors in a safe and consistent way.
 */
class SafeShapeExtractor {
public:
    /**
     *  Construct the extractor, creating a flag alias that indicates failure in the input centroid
     *  by linking to the slot shape flag.
     *
     *  @param[out] schema   Schema to which the alias should be added.  The "slot_Shape" alias
     *                       must already be present in the Schema's AliasMap.
     *  @param[in]  name     The name of the algorithm; the flag alias added will be
     *                       "<name>_flag_badShape".
     */
    SafeShapeExtractor(afw::table::Schema& schema, std::string const& name);

    /**
     *  Extract a shape from the given record.
     *
     *  We use the Shape slot if it is not NaN, and throw MeasurementError if it is (with only the general
     *  failure code attached to the MeasurementError).
     *
     *  If the Shape slot is not defined, we throw FatalAlgorithmError, as this indicates a configuration
     *  problem.
     *
     *  If the Shape slot value is NaN and is *not* flagged (or there is no Shape slot flag), we throw
     *  RuntimeError, which should cause the measurement framework to log a warning and set the current
     *  algorithm's general failure flag if it is allowed to propagate out of the algorithm implementation.
     *
     *  If the Shape slot is flagged and we nevertheless obtain a usable ellipse, we set the current
     *  algorithm's general failure flag, but return the ellipse as well, allowing it to continue while
     *  indicating that the result may not be reliable.  A singular ellipse (i.e. one with a non-positive
     *  quadrupole matrix determinant) is treated the same as a NaN ellipse; it is not considered usable.
     */
    afw::geom::ellipses::Quadrupole operator()(afw::table::SourceRecord& record,
                                               FlagHandler const& flags) const;

private:
    std::string _name;
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_InputUtilities_h_INCLUDED
