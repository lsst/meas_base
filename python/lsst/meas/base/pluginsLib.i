
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

// Centroid algorithms
%{
#include "lsst/meas/base/SdssCentroid.h"
%}

%include "lsst/meas/base/SdssCentroid.h"

%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SdssCentroidAlgorithm,
                           SdssCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/GaussianCentroid.h"
%template(apply) lsst::meas::base::GaussianCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::GaussianCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, GaussianCentroidAlgorithm, GaussianCentroidControl,
                           FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/NaiveCentroid.h"
%template(apply) lsst::meas::base::NaiveCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::NaiveCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, NaiveCentroidAlgorithm,
                           NaiveCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/SdssCentroid.h"
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SdssCentroidAlgorithm,
                           SdssCentroidControl, FootprintCentroidInput, CentroidComponent)

// Flux algorithms
%include "lsst/meas/base/PsfFlux.h"
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, PsfFluxAlgorithm, PsfFluxControl, FootprintCentroidInput,
                           FluxComponent)

%include "lsst/meas/base/SincFlux.h"
%template(apply) lsst::meas::base::SincFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SincFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SincFluxAlgorithm, SincFluxControl,
                           FootprintCentroidInput, FluxComponent)

%include "lsst/meas/base/GaussianFlux.h"
%template(apply) lsst::meas::base::GaussianFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::GaussianFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, GaussianFluxAlgorithm, GaussianFluxControl,
                           FootprintCentroidShapeInput, FluxComponent)

%include "lsst/meas/base/NaiveFlux.h"
%template(apply) lsst::meas::base::NaiveFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::NaiveFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, NaiveFluxAlgorithm, NaiveFluxControl,
                           FootprintCentroidInput, FluxComponent)

%include "lsst/meas/base/PeakLikelihoodFlux.h"
%template(apply) lsst::meas::base::PeakLikelihoodFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PeakLikelihoodFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, PeakLikelihoodFluxAlgorithm, PeakLikelihoodFluxControl,
                           FootprintCentroidInput, FluxComponent)

// Shape algorithms

%include "lsst/meas/base/SdssShape.h"
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<double>;
%wrapMeasurementAlgorithm4(lsst::meas::base, SdssShapeAlgorithm, SdssShapeControl, FootprintCentroidInput,
                          ShapeComponent, CentroidComponent, FluxComponent, SdssShapeExtras)

// Miscellaneous algorithms

%include "lsst/meas/base/PixelFlags.h"
%template(apply) lsst::meas::base::PixelFlagsAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PixelFlagsAlgorithm::apply<double>;
%wrapMeasurementAlgorithm0(lsst::meas::base, PixelFlagsAlgorithm, PixelFlagsControl, FootprintCentroidInput)
%include "lsst/meas/base/PixelFlags.h"

%include "ApertureFlux.i"
