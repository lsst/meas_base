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
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/table.h"
#include "lsst/pex/exceptions.h"

namespace lsst {
namespace meas {
namespace base {

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
 *  are contiguous in memory. It is good practice to ensure that they are
 *  equal in size: this may be conveniently achieved by calling
 *  `BaseTransform::checkCatalogSize()`.
 *
 *  `operator()` may throw `LengthError` if the transformation is impossible
 *  to complete. In this case, the contents of `outputCatalog` is not
 *  guaranteed.
 */
class BaseTransform {
public:
    explicit BaseTransform(std::string const& name) : _name(name) {}
    virtual ~BaseTransform() {}
    virtual void operator()(afw::table::SourceCatalog const& inputCatalog,
                            afw::table::BaseCatalog& outputCatalog, afw::geom::SkyWcs const& wcs,
                            afw::image::PhotoCalib const& photoCalib) const = 0;

protected:
    /**
     *  @brief Ensure that catalogs have the same size.
     *
     *  @param[in]  cat1         Catalog for comparison
     *  @param[in]  cat2         Catalog for comparison
     *  @throws     LengthError  Catalog sizes do not match
     */
    void checkCatalogSize(afw::table::BaseCatalog const& cat1, afw::table::BaseCatalog const& cat2) const {
        if (cat1.size() != cat2.size()) {
            throw LSST_EXCEPT(pex::exceptions::LengthError, "Catalog size mismatch");
        }
    }
    std::string _name;
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_Transform_h_INCLUDED
