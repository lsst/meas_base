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

%{
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/afw/table.h"
%}

%import "lsst/afw/table/tableLib.i"

%immutable lsst::meas::base::FlagDefinition::name;
%immutable lsst::meas::base::FlagDefinition::doc;

%declareNumPyConverters(lsst::meas::base::CentroidCov);
%declareNumPyConverters(lsst::meas::base::ShapeCov);
%declareNumPyConverters(lsst::meas::base::ShapeTrMatrix);

%declareFunctorKey(FluxResult, lsst::meas::base::FluxResult)
%shared_ptr(lsst::meas::base::FluxResultKey)

%declareFunctorKey(MagResult, lsst::meas::base::MagResult)
%shared_ptr(lsst::meas::base::MagResultKey)

%declareFunctorKey(CentroidResult, lsst::meas::base::CentroidResult)
%shared_ptr(lsst::meas::base::CentroidResultKey)

%declareFunctorKey(ShapeResult, lsst::meas::base::ShapeResult)
%shared_ptr(lsst::meas::base::ShapeResultKey)

%include "lsst/meas/base/FluxUtilities.h"
%include "lsst/meas/base/CentroidUtilities.h"
%include "lsst/meas/base/ShapeUtilities.h"
%include "lsst/meas/base/FlagHandler.h"
%include "lsst/meas/base/InputUtilities.h"
