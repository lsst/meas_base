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
%include "lsst/meas/base/exceptions.h"
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

%instantiateInput(FootprintInput)
%instantiateInput(FootprintCentroidInput)
%instantiateInput(FootprintCentroidShapeInput)

%define %instantiateResult0(NAMESPACE, ALGORITHM)
%template(ALGORITHM##FlagsComponent) lsst::meas::base::FlagsComponent<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##FlagsComponentMapper) lsst::meas::base::FlagsComponentMapper<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##Result0) lsst::meas::base::Result0<
    NAMESPACE::ALGORITHM
    >;
%template(ALGORITHM##ResultMapper0) lsst::meas::base::ResultMapper0<
    NAMESPACE::ALGORITHM
    >;
%enddef

%define %instantiateResult1(NAMESPACE, ALGORITHM, T1)
%template(ALGORITHM##FlagsComponent) lsst::meas::base::FlagsComponent<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##FlagsComponentMapper) lsst::meas::base::FlagsComponentMapper<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##Result1) lsst::meas::base::Result1<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1
    >;
%template(ALGORITHM##ResultMapper1) lsst::meas::base::ResultMapper1<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Mapper
    >;
%enddef

%define %instantiateResult2(NAMESPACE, ALGORITHM, T1, T2)
%template(ALGORITHM##FlagsComponent) lsst::meas::base::FlagsComponent<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##FlagsComponentMapper) lsst::meas::base::FlagsComponentMapper<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##Result2) lsst::meas::base::Result2<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1,
    lsst::meas::base::T2
    >;
%template(ALGORITHM##ResultMapper2) lsst::meas::base::ResultMapper2<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Mapper,
    lsst::meas::base::T2##Mapper
    >;
%enddef

%define %instantiateResult3(NAMESPACE, ALGORITHM, T1, T2, T3)
%template(ALGORITHM##FlagsComponent) lsst::meas::base::FlagsComponent<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##FlagsComponentMapper) lsst::meas::base::FlagsComponentMapper<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##Result3) lsst::meas::base::Result3<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1,
    lsst::meas::base::T2,
    lsst::meas::base::T3
    >;
%template(ALGORITHM##ResultMapper3) lsst::meas::base::ResultMapper3<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Mapper,
    lsst::meas::base::T2##Mapper,
    lsst::meas::base::T3##Mapper
    >;
%enddef

%define %instantiateResult4(NAMESPACE, ALGORITHM, T1, T2, T3, T4)
%template(ALGORITHM##FlagsComponent) lsst::meas::base::FlagsComponent<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##FlagsComponentMapper) lsst::meas::base::FlagsComponentMapper<NAMESPACE::ALGORITHM>;
%template(ALGORITHM##Result4) lsst::meas::base::Result4<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1,
    lsst::meas::base::T2,
    lsst::meas::base::T3,
    lsst::meas::base::T4
    >;
%template(ALGORITHM##ResultMapper4) lsst::meas::base::ResultMapper4<
    NAMESPACE::ALGORITHM,
    lsst::meas::base::T1##Mapper,
    lsst::meas::base::T2##Mapper,
    lsst::meas::base::T3##Mapper,
    lsst::meas::base::T4##Mapper
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

%define %wrapMeasurementAlgorithm0(NAMESPACE, ALGORITHM, CONTROL, INPUT)
%instantiateResult0(NAMESPACE, ALGORITHM)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##Result0, ALGORITHM##ResultMapper0)
%enddef

%define %wrapMeasurementAlgorithm1(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1)
%instantiateResult1(NAMESPACE, ALGORITHM, T1)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##Result1, ALGORITHM##ResultMapper1)
%enddef

%define %wrapMeasurementAlgorithm2(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1, T2)
%instantiateResult2(NAMESPACE, ALGORITHM, T1, T2)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##Result2, ALGORITHM##ResultMapper2)
%enddef

%define %wrapMeasurementAlgorithm3(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1, T2, T3)
%instantiateResult3(NAMESPACE, ALGORITHM, T1, T2, T3)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##Result3, ALGORITHM##ResultMapper3)
%enddef

%define %wrapMeasurementAlgorithm4(NAMESPACE, ALGORITHM, CONTROL, INPUT, T1, T2, T3, T4)
%instantiateResult4(NAMESPACE, ALGORITHM, T1, T2, T3, T4)
%wrapMeasurementAlgorithmEx(NAMESPACE, ALGORITHM, CONTROL, INPUT,
                            ALGORITHM##Result4, ALGORITHM##ResultMapper4)
%enddef

%include "lsst/meas/base/PsfFlux.h"
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PsfFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, PsfFluxAlgorithm, PsfFluxControl, FootprintCentroidInput, FluxComponent)

%include "lsst/meas/base/SdssShape.h"
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssShapeAlgorithm::apply<double>;
%wrapMeasurementAlgorithm4(lsst::meas::base, SdssShapeAlgorithm, SdssShapeControl, FootprintCentroidInput,
                          ShapeComponent, CentroidComponent, FluxComponent, SdssShapeExtras)

%include "lsst/meas/base/SincFlux.h"
%template(apply) lsst::meas::base::SincFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SincFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SincFluxAlgorithm, SincFluxControl, FootprintCentroidInput, FluxComponent)

%include "lsst/meas/base/GaussianFlux.h"
%template(apply) lsst::meas::base::GaussianFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::GaussianFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, GaussianFluxAlgorithm, GaussianFluxControl, FootprintCentroidInput, FluxComponent)

%include "lsst/meas/base/GaussianCentroid.h"
%template(apply) lsst::meas::base::GaussianCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::GaussianCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, GaussianCentroidAlgorithm, GaussianCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/NaiveFlux.h"
%template(apply) lsst::meas::base::NaiveFluxAlgorithm::apply<float>;
%template(apply) lsst::meas::base::NaiveFluxAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, NaiveFluxAlgorithm, NaiveFluxControl, FootprintCentroidInput, FluxComponent)

%include "lsst/meas/base/NaiveCentroid.h"
%template(apply) lsst::meas::base::NaiveCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::NaiveCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, NaiveCentroidAlgorithm, NaiveCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/SdssCentroid.h"
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<float>;
%template(apply) lsst::meas::base::SdssCentroidAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, SdssCentroidAlgorithm, SdssCentroidControl, FootprintCentroidInput, CentroidComponent)

%include "lsst/meas/base/PixelFlags.h"
%template(apply) lsst::meas::base::PixelFlagsAlgorithm::apply<float>;
%template(apply) lsst::meas::base::PixelFlagsAlgorithm::apply<double>;
%wrapMeasurementAlgorithm0(lsst::meas::base, PixelFlagsAlgorithm, PixelFlagsControl, FootprintCentroidInput)
%include "lsst/meas/base/PixelFlags.h"

%include "lsst/meas/base/Classification.h"
%template(apply) lsst::meas::base::ClassificationAlgorithm::apply<float>;
%template(apply) lsst::meas::base::ClassificationAlgorithm::apply<double>;
%wrapMeasurementAlgorithm1(lsst::meas::base, ClassificationAlgorithm, ClassificationControl, ClassificationInput, ClassificationExtras)
%include "lsst/meas/base/Classification.h"
