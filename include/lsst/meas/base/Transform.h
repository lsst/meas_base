// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST.
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

#ifndef LSST_MEAS_BASE_Transform_h_INCLUDED
#define LSST_MEAS_BASE_Transform_h_INCLUDED

/**
 *  @file lsst/meas/base/Transform.h
 *  This defines the base of measurement transformations.
 */

#include <string>
#include "lsst/afw/image.h"
#include "lsst/afw/table.h"

namespace lsst { namespace meas { namespace base {

/**
 *  Abstract base class for all C++ measurement transformations
 *
 *  Measurement plugins return results in raw, uncalibrated units (eg fluxes
 *  or positions in pixels). The transformation system provides a mechanism
 *  for post-processing those results into a calibrated form (magnitudes,
 *  celestial coordinates, etc).
 *
 *  A measurement transformation should derive from BaseTransform. It should
 *  implement a constructor which takes three arguments:
 *
 *  - A `Control` object describing the configuration of the measurement
 *    plugin.
 *  - The name of the measurement plugin whose outputs are to be transformed
 *    (`std::string`);
 *  - An `lsst::afw::table::SchemaMapper` which links the input and output
 *    catalogs;
 *
 *  The constructor should use the SchemaMapper to map fields from the input
 *  to output schemas and add additional keys to the output as required. For
 *  example:
 *
 *  @dontinclude SillyCentroid.h
 *  @skip SillyTransform
 *  @until }
 *
 *  Derived classes should also implement `operator()` following the interface
 *  below. This will be called with a catalog containing the results of the
 *  measurement plugin and a catalog to be populated with transformed
 *  quantities, as well as WCS and calibration information. For example:
 *
 *  @skip operator()
 *  @until // operator()
 *
 *  Note that it is safe to assume that both catalogs passed to `operator()`
 *  are contiguous in memory.
 *
 */
class BaseTransform {
public:

    explicit BaseTransform(std::string const & name) : _name(name) {}
    virtual ~BaseTransform() {}
    virtual void operator()(afw::table::SourceCatalog const & inputCatalog,
                            afw::table::BaseCatalog & outputCatalog,
                            afw::image::Wcs const & wcs,
                            afw::image::Calib const & calib) const = 0;

protected:
    std::string _name;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_Transform_h_INCLUDED
