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
%import "lsst/afw/image/Exposure.h"

%immutable lsst::meas::base::FlagDef::name;
%immutable lsst::meas::base::FlagDef::doc;

%include "lsst/meas/base/Results.h"
%include "lsst/meas/base/ResultMappers.h"
%include "lsst/meas/base/Inputs.h"

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

%pythoncode %{
    FlagsResult = {}
    FlagsResultMapper = {}
%}

%define %instantiateFlags(N)
%template(FlagsResult ## N) lsst::meas::base::FlagsResult<N>;
%template(FlagsResultMapper ## N) lsst::meas::base::FlagsResultMapper<N>;
%pythoncode %{
    FlagsResult[N] = FlagsResult ## N
    FlagsResultMapper[N] = FlagsResultMapper ## N
%}
%enddef

%instantiateFlags(0)
%instantiateFlags(1)
%instantiateFlags(2)
%instantiateFlags(3)
%instantiateFlags(4)
%instantiateFlags(5)
%instantiateFlags(6)
%instantiateFlags(7)
%instantiateFlags(8)

%define %instantiateSimpleResult1(NAMESPACE, ALGORITHM, T1)
%template(ALGORITHM##SimpleResult1) lsst::meas::base::SimpleResult1<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Result
    >;
%template(ALGORITHM##SimpleResultMapper1) lsst::meas::base::SimpleResult1<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##ResultMapper
    >;
%enddef

%define %instantiateSimpleResult2(NAMESPACE, ALGORITHM, T1, T2)
%template(ALGORITHM##SimpleResult2) lsst::meas::base::SimpleResult2<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Result,
    lsst::meas::base::T2##Result
    >;
%template(ALGORITHM##SimpleResultMapper2) lsst::meas::base::SimpleResult2<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##ResultMapper,
    lsst::meas::base::T2##ResultMapper
    >;
%enddef

%define %instantiateSimpleResult3(NAMESPACE, ALGORITHM, T1, T2, T3)
%template(ALGORITHM##SimpleResult3) lsst::meas::base::SimpleResult3<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Result,
    lsst::meas::base::T2##Result,
    lsst::meas::base::T3##Result
    >;
%template(ALGORITHM##SimpleResultMapper3) lsst::meas::base::SimpleResult3<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##ResultMapper,
    lsst::meas::base::T2##ResultMapper,
    lsst::meas::base::T3##ResultMapper
    >;
%enddef

%define %wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT, RESULT, RESULT_MAPPER)
// Turn C++ typedefs into equivalent Python attributes - it's a shame Swig doesn't do this for us,
// or even give us a way to look up the Python class name for a C++ class we've already wrapped.
%pythoncode %{
ALGORITHM.Control = CONTROL
ALGORITHM.Input = INPUT
ALGORITHM.Result = RESULT
ALGORITHM.ResultMapper = RESULT_MAPPER
%}
%enddef

%define %wrapMeasurementAlgorithm1(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1)
%instantiateSimpleResult1(NAMESPACE, ALGORITHM, T1)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##SimpleResult1, ALGORITHM##SimpleResultMapper1)
%enddef

%define %wrapMeasurementAlgorithm2(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1, T2)
%instantiateSimpleResult2(NAMESPACE, ALGORITHM, T1, T2)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##SimpleResult2, ALGORITHM##SimpleResultMapper2)
%enddef

%define %wrapMeasurementAlgorithm3(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1, T2, T3)
%instantiateSimpleResult3(NAMESPACE, ALGORITHM, T1, T2, T3)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##SimpleResult3, ALGORITHM##SimpleResultMapper3)
%enddef

%include "lsst/meas/base/PsfFlux.h"
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<double>;
%template(applyN) lsst::meas::base::PsfFluxAlgorithm::applyN<float>;
%template(applyN) lsst::meas::base::PsfFluxAlgorithm::applyN<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, PsfFluxAlgorithm, NullControl, AlgorithmInput2, Flux)

%include "lsst/meas/base/SdssShape.h"
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<double>;
%wrapMeasurementAlgorithm3(lsst::meas::base, SdssShapeAlgorithm, SdssShapeControl, AlgorithmInput2,
                           Shape, Centroid, Flux)
