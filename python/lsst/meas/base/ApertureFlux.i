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

%{
#include "lsst/meas/base/ApertureFlux.h"
#include "lsst/meas/base/CircularApertureFlux.h"
%}

%feature("abstract") lsst::meas::base::ApertureFluxAlgorithm;
%include "lsst/meas/base/ApertureFlux.h"

%pythoncode %{
ApertureFluxAlgorithm.Result = ApertureFluxResult;
ApertureFluxAlgorithm.Control = ApertureFluxControl;
%}
%template(computeNaiveFlux) lsst::meas::base::ApertureFluxAlgorithm::computeNaiveFlux<float>;
%template(computeNaiveFlux) lsst::meas::base::ApertureFluxAlgorithm::computeNaiveFlux<double>;
%template(computeSincFlux) lsst::meas::base::ApertureFluxAlgorithm::computeSincFlux<float>;
%template(computeSincFlux) lsst::meas::base::ApertureFluxAlgorithm::computeSincFlux<double>;
%template(computeFlux) lsst::meas::base::ApertureFluxAlgorithm::computeFlux<float>;
%template(computeFlux) lsst::meas::base::ApertureFluxAlgorithm::computeFlux<double>;

%feature("notabstract") lsst::meas::base::CircularApertureFluxAlgorithm;
%include "lsst/meas/base/CircularApertureFlux.h"
