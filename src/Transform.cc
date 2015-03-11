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

#include "lsst/afw/table.h"
#include "lsst/meas/base/Transform.h"

namespace lsst { namespace meas { namespace base {

FluxTransform::FluxTransform(
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    BaseTransform{name}
{
    // Map the flag through to the output
    mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(name + "_flag").key);

    // Add keys for the magnitude and associated error
    _magKey = mapper.editOutputSchema().addField<double>(name + "_mag", "Magnitude");
    _magErrKey = mapper.editOutputSchema().addField<double>(name + "_magErr", "Magnitude Error");
}

void FluxTransform::operator()(
    afw::table::SourceCatalog const & inputCatalog,
    afw::table::BaseCatalog & outputCatalog,
    afw::image::Wcs const & wcs,
    afw::image::Calib const & calib
) const {
    afw::table::Key<double> fluxKey = inputCatalog.getSchema().find<double>(_name + "_flux").key;
    afw::table::Key<double> fluxSigmaKey = inputCatalog.getSchema().find<double>(_name + "_fluxSigma").key;

    afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
    afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();
    double mag, magErr;
    for (; inSrc < inputCatalog.end() && outSrc < outputCatalog.end(); ++inSrc, ++outSrc) {
        boost::tie(mag, magErr) = calib.getMagnitude(inSrc->get(fluxKey), inSrc->get(fluxSigmaKey));
        outSrc->set(_magKey, mag);
        outSrc->set(_magErrKey, magErr);
    }
}

}}} // namespace lsst::meas::base
