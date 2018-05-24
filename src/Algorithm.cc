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

#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/Algorithm.h"

namespace lsst {
namespace meas {
namespace base {

void SingleFrameAlgorithm::measureN(afw::table::SourceCatalog const& measCat,
                                    afw::image::Exposure<float> const& exposure) const {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "measureN not implemented for this algorithm");
}

void ForcedAlgorithm::measureNForced(afw::table::SourceCatalog const& measCat,
                                     afw::image::Exposure<float> const& exposure,
                                     afw::table::SourceCatalog const& refRecord,
                                     afw::geom::SkyWcs const& refWcs) const {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "measureN not implemented for this algorithm");
}

}  // namespace base
}  // namespace meas
}  // namespace lsst
