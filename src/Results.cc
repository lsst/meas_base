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

FluxComponent::FluxComponent() :
    flux(std::numeric_limits<Flux>::quiet_NaN()),
    fluxSigma(std::numeric_limits<ErrElement>::quiet_NaN())
{}

CentroidComponent::CentroidComponent() :
    x(std::numeric_limits<CentroidElement>::quiet_NaN()),
    y(std::numeric_limits<CentroidElement>::quiet_NaN()),
    xSigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    ySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    x_y_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}

ShapeComponent::ShapeComponent() :
    xx(std::numeric_limits<ShapeElement>::quiet_NaN()),
    yy(std::numeric_limits<ShapeElement>::quiet_NaN()),
    xy(std::numeric_limits<ShapeElement>::quiet_NaN()),
    xxSigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    yySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    xySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    xx_yy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    xx_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN()),
    yy_xy_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}

}}} // namespace lsst::meas::base
