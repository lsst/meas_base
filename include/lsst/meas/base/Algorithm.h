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

#ifndef LSST_MEAS_BASE_Algorithm_h_INCLUDED
#define LSST_MEAS_BASE_Algorithm_h_INCLUDED

#include "lsst/afw/table/fwd.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/base/FlagHandler.h"

namespace lsst { namespace meas { namespace base {

/**
 *  Ultimate abstract base class for all C++ measurement algorithms
 *
 *  New algorithms should not inherit directly from this class.
 */
class BaseAlgorithm {
public:

    /**
     *  Handle an exception thrown by the current algorithm by setting flags in the given
     *  record.
     *
     *  fail() is called by the measurement framework when an exception is allowed to
     *  propagate out of one the algorithm's measure() methods.  It should generally set
     *  both a general failure flag for the algorithm as well as a specific flag
     *  indicating the error condition, if possible.  To aid in this, if the exception was
     *  an instance of MeasurementError, it will be passed in, carrying information about
     *  what flag to set.
     *
     *  An algorithm can also to chose to set flags within its own measure() methods, and
     *  then just return, rather than throw an exception.  However, fail() should be
     *  implemented even when all known failure modes do not throw exceptions, to ensure
     *  that unexpected exceptions thrown in lower-level code are properly handled.
     */
    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const = 0;

    virtual ~BaseAlgorithm() {}

};

/**
 *  Base class for algorithms that measure the properties of sources on single image.
 *
 *  SingleFrameAlgorithm defines the interface used in measuring both on single exposure images
 *  and on coadds.
 *
 *  In addition to the virtual methods defined here, SingleFrameAlgorithm also puts requirements
 *  on constructor signatures; see the wrapSingleFrameAlgorithm Python function for more
 *  information.
 */
class SingleFrameAlgorithm : public virtual BaseAlgorithm {
public:

    /**
     *  Called to measure a single child source in an image.
     *
     *  Before this method is called, all neighbors will be replaced with noise, using the
     *  outputs of the deblender.  Outputs should be saved in the given SourceRecord,
     *  which can also be used to obtain centroid (see SafeCentroidExtractor) and shape
     *  (see SafeShapeExtractor) information.
     */
    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const = 0;

    /**
     *  Called to simultaneously measure all children in a deblend family, in a single image.
     *
     *  Outputs should be saved in the given SourceCatalog, which can also be used to
     *  obtain centroid (see SafeCentroidExtractor) and shape (see SafeShapeExtractor)
     *  information.
     *
     *  The default implementation simply throws an exception, indicating that simultaneous
     *  measurement is not supported.
     */
    virtual void measureN(
        afw::table::SourceCatalog const & measCat,
        afw::image::Exposure<float> const & exposure
    ) const;

};

/**
 *  Base class for algorithms that measure the properties of sources on one image, using
 *  previous measurements on another image to hold certain quantities fixed.
 *
 *  ForcedAlgorithm can be used when measuring both on single exposures and coadds, and is typically
 *  used to measure colors by holding centroids, apertures, or other model parameters fixed while
 *  allowing amplitudes to vary.
 *
 *  In addition to the virtual methods defined here, ForcedAlgorithm also puts requirements
 *  on constructor signatures; see the wrapForcedAlgorithm Python function for more
 *  information.
 *
 *  Most algorithms will not need to make use of the reference record or WCS, as transformed
 *  centroids and shapes should generally be available via the slots in measRecord, but these
 *  are made available for algorithms that need to transform more complex information.  If that
 *  is the case, the algorithm may want to inherit from SimpleAlgorithm instead of inheriting
 *  from ForcedAlgorithm directly.
 */
class ForcedAlgorithm : public virtual BaseAlgorithm {
public:

    /**
     *  Called to measure a single child source in an image.
     *
     *  Before this method is called, all neighbors will be replaced with noise, using the
     *  outputs of the deblender.  Outputs should be saved in the given SourceRecord,
     *  which can also be used to obtain centroid (see SafeCentroidExtractor) and shape
     *  (see SafeShapeExtractor) information.
     */
    virtual void measureForced(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceRecord const & refRecord,
        afw::image::Wcs const & refWcs
    ) const = 0;

    /**
     *  Called to simultaneously measure all children in a deblend family, in a single image.
     *
     *  Outputs should be saved in the given SourceCatalog, which can also be used to
     *  obtain centroid (see SafeCentroidExtractor) and shape (see SafeShapeExtractor)
     *  information.
     *
     *  The default implementation simply throws an exception, indicating that simultaneous
     *  measurement is not supported.
     */
    virtual void measureNForced(
        afw::table::SourceCatalog const & measCat,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceCatalog const & refRecord,
        afw::image::Wcs const & refWcs
    ) const;

};

/**
 *  An abstract base classes for which the same implementation can be used for both SingleFrameAlgorithm
 *  and ForcedAlgorithm.
 *
 *  SimpleAlgorithm allows a ForcedAlgorithm to be defined using the measure() and measureN() signatures
 *  of SingleFrameAlgorithm.  It should be used for any algorithm for whichi the forced version of the
 *  algorithm does not require anything to be transformed beyond the centroid and shape, and for which
 *  the the parameters being fit are the same for both single-frame and forced measurement.  This should
 *  be the case for all flux algorithms that don't involve fitting any additional model parameters.  It
 *  can also be used for centroid and shape algorithms, where having a version that can re-measure values in
 *  forced mode may be useful for diagnostic purposes even if it is not useful for science.
 */
class SimpleAlgorithm : public SingleFrameAlgorithm, public ForcedAlgorithm {
public:

    virtual void measureForced(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceRecord const & refRecord,
        afw::image::Wcs const & refWcs
    ) const {
        measure(measRecord, exposure);
    }

    virtual void measureNForced(
        afw::table::SourceCatalog const & measCat,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceCatalog const & refRecord,
        afw::image::Wcs const & refWcs
    ) const {
        measureN(measCat, exposure);
    }

};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_Algorithm_h_INCLUDED
