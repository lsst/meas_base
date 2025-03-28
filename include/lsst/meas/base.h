// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
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

#ifndef LSST_MEAS_base_h_INCLUDED
#define LSST_MEAS_base_h_INCLUDED

#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/PsfFlux.h"
#include "lsst/meas/base/SdssCentroid.h"
#include "lsst/meas/base/SdssShape.h"
#include "lsst/meas/base/PixelFlags.h"
#include "lsst/meas/base/GaussianFlux.h"
#include "lsst/meas/base/PeakLikelihoodFlux.h"
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/ScaledApertureFlux.h"
#include "lsst/meas/base/CircularApertureFlux.h"
#include "lsst/meas/base/Blendedness.h"

// These are necessary to build Swig modules that %import meas/base/baseLib.i,
// so it's neighborly to include them here so downstream code can just
// #include "lsst/meas/base.h" in their Swig wrappers, instead of guessing
// what else they might need.
#include "lsst/afw/detection.h"
#include "lsst/afw/math.h"

#endif  // !LSST_MEAS_base_h_INCLUDED
