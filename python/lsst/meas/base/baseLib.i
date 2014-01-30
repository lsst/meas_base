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

%define baseLib_DOCSTRING
"
Basic routines to talk to lsst::meas::base classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.base", docstring=baseLib_DOCSTRING) baseLib

%{
#include "lsst/pex/logging.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/math.h"
#include "lsst/afw/table.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/meas/base.h"
#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_BASE_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%include "std_complex.i"

%include "ndarray.i"

%declareNumPyConverters(lsst::meas::base::CentroidCov);

%lsst_exceptions();

%include "std_vector.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/pex/config.h"

%include "lsst/meas/base/Results.h"
%include "lsst/meas/base/ResultMappers.h"
%include "lsst/afw/image/Exposure.h"
%include "lsst/meas/base/Inputs.h"
%include "lsst/meas/base/PsfFlux.h"
%include "lsst/meas/base/SdssShape.h"

%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<double>;
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<double>;

%define %instantiateInput(T)
%ignore std::vector<lsst::meas::base::T>::vector(size_type);
%ignore std::vector<lsst::meas::base::T>::resize(size_type);
%template(T ## Vector) std::vector<lsst::meas::base::T>;
%pythoncode %{
T.Vector = T ## Vector
%}
%enddef

%instantiateInput(AlgorithmInput1)
%instantiateInput(AlgorithmInput2)
%instantiateInput(AlgorithmInput3)

// Turn C++ typedefs into equivalent Python attributes - it's a shame Swig doesn't do this for us.
%pythoncode %{

PsfFluxAlgorithm.Control = NullControl
PsfFluxAlgorithm.Result = FluxResult
PsfFluxAlgorithm.ResultMapper = FluxResultMapper
PsfFluxAlgorithm.Input = AlgorithmInput2

SdssShapeAlgorithm.Control = SdssShapeControl
SdssShapeAlgorithm.Result = SdssShapeResult
SdssShapeAlgorithm.ResultMapper = SdssShapeResultMapper
SdssShapeAlgorithm.Input = AlgorithmInput2

%}
