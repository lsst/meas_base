// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * See COPYRIGHT file at the top of the source tree.
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

#ifndef LSST_MEAS_BASE_PYTHON_H
#define LSST_MEAS_BASE_PYTHON_H

#include "pybind11/pybind11.h"

#include <string>

#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/table/fwd.h"
#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/SchemaMapper.h"

namespace lsst {
namespace meas {
namespace base {
namespace python {

namespace py = pybind11;
using namespace pybind11::literals;

/**
 * Wrap the standard algorithm constructor.
 *
 * @tparam Algorithm The algorithm class.
 * @tparam PyAlg The `pybind11::class_` class corresponding to `Algorithm`.
 *
 * @param[in,out] cls The pybind11 wrapper for `Algorithm`.
 */
template <class Algorithm, class PyAlg>
typename std::enable_if<!std::is_abstract<Algorithm>::value, void>::type declareAlgorithmConstructor(PyAlg & cls) {
    cls.def(py::init<typename Algorithm::Control const &,
                     std::string const &,
                     afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
}
/**
 * Dummy function for not wrapping the constructor of an abstract base class.
 *
 * Pybind11 cannot wrap such constructors, and there is no reason to call them
 * from Python anyway.
 *
 * @tparam Algorithm The algorithm class.
 * @tparam PyAlg The `pybind11::class_` class corresponding to `Algorithm`.
 *
 * @param[in,out] cls The pybind11 wrapper for `Algorithm`.
 */
template <class Algorithm, class PyAlg>
typename std::enable_if<std::is_abstract<Algorithm>::value, void>::type declareAlgorithmConstructor(PyAlg & cls) {
}

/**
 * Wrap the implicit API used by meas_base's algorithms.
 *
 * This function only initializes constructors, fields, and methods common to
 * all Algorithms.
 *
 * @tparam Algorithm The algorithm class.
 * @tparam PyAlg The `pybind11::class_` class corresponding to `Algorithm`.
 *
 * @param[in,out] clsAlgorithm The pybind11 wrapper for `Algorithm`.
 */
template <class Algorithm, class PyAlg>
void declareAlgorithm(PyAlg & clsAlgorithm) {
    /* Member types and enums */

    /* Constructors */
    declareAlgorithmConstructor<Algorithm>(clsAlgorithm);

    /* Operators */

    /* Members */
    clsAlgorithm.def("fail", &Algorithm::fail, "measRecord"_a, "error"_a=NULL);
    clsAlgorithm.def("measure", &Algorithm::measure, "record"_a, "exposure"_a);
}

/**
 * Wrap the implicit API used by meas_base's algorithm-control pairs (no transform).
 *
 * This function only initializes constructors, fields, and methods common to
 * all Algorithms and Controls.
 *
 * @tparam Algorithm The algorithm class.
 * @tparam Control The control class. Must equal `Algorithm::Control` and
 *                 `Transform::Control`.
 * @tparam PyAlg The `pybind11::class_` class corresponding to `Algorithm`.
 * @tparam PyCtrl The `pybind11::class_` class corresponding to `Control`.
 *
 * @param[in,out] clsAlgorithm,clsControl The pybind11 wrappers
 *                             for the respective C++ classes.
 */
template <class Algorithm, class Control, class PyAlg, class PyCtrl>
void declareAlgorithm(PyAlg & clsAlgorithm, PyCtrl & clsControl) {
    declareAlgorithm<Algorithm>(clsAlgorithm);

    /* Member types and enums */

    /* Constructors */
    clsControl.def(py::init<>());

    /* Operators */

    /* Members */
}

/**
 * Wrap the implicit API used by meas_base's algorithm-control-transform
 * triads.
 *
 * This function only initializes constructors, fields, and methods common to
 * all Algorithms, Controls, and Transforms.
 *
 * @tparam Algorithm The algorithm class.
 * @tparam Control The control class. Must equal `Algorithm::Control` and
 *                 `Transform::Control`.
 * @tparam Transform The transform class.
 * @tparam PyAlg The `pybind11::class_` class corresponding to `Algorithm`.
 * @tparam PyCtrl The `pybind11::class_` class corresponding to `Control`.
 * @tparam PyXform The `pybind11::class_` class corresponding to `Transform`.
 *
 * @param[in,out] clsAlgorithm,clsControl,clsTransform The pybind11 wrappers
 *                             for the respective C++ classes.
 */
template <class Algorithm, class Control, class Transform, class PyAlg, class PyCtrl, class PyXform>
void declareAlgorithm(PyAlg & clsAlgorithm, PyCtrl & clsControl, PyXform & clsTransform) {
    declareAlgorithm<Algorithm, Control>(clsAlgorithm, clsControl);

    /* Member types and enums */

    /* Constructors */
    clsTransform.def(py::init<typename Transform::Control const &,
                              std::string const &,
                              afw::table::SchemaMapper &>(),
                     "ctrl"_a, "name"_a, "mapper"_a);

    /* Operators */
    clsTransform.def("__call__", [](Transform const & self,
                                    afw::table::SourceCatalog const & inputCatalog,
                                    afw::table::BaseCatalog & outputCatalog,
                                    afw::image::Wcs const & wcs,
                                    afw::image::Calib const &calib) {
            return self(inputCatalog, outputCatalog, wcs, calib);
        }, "inputCatalog"_a, "outputCatalog"_a, "wcs"_a, "calib"_a);

    /* Members */
}

}}}}     // lsst::meas::base::

#endif

