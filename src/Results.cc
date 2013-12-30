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

#include "lsst/meas/base/Results.h"

namespace lsst { namespace meas { namespace base {

FluxAlgorithmResult::FluxAlgorithmResult() :
    value(std::numeric_limits<Flux>::quiet_NaN()),
    err(std::numeric_limits<FluxErr>::quiet_NaN())
{}

CentroidAlgorithmResult::CentroidAlgorithmResult() :
    value(std::numeric_limits<double>::quiet_NaN()),
    cov(CentroidCov::Constant(std::numeric_limits<float>::quiet_NaN()))
{}

ShapeAlgorithmResult::ShapeAlgorithmResult() :
    value(std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN()),
    cov(ShapeCov::Constant(std::numeric_limits<float>::quiet_NaN()))
{}

}}} // namespace lsst::meas::base
