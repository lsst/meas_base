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

#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

BaseAlgorithmMapper::BaseAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix
) :
    _flag(
        schema.addField(
            afw::table::Field<afw::table::Flag>(
                prefix + ".flag",
                "general failure flag for " + prefix + " measurement"
            ),
            true // replace existing fields if present
        )
    )
{}

void BaseAlgorithmMapper::fail(afw::table::BaseRecord & record) {
    record.set(_flag, true);
}

FluxAlgorithmMapper::FluxAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix
) : BaseAlgorithmMapper(schema, prefix),
    _value(schema.addField(
               afw::table::Field<Flux>(prefix + ".value", "measured flux", "dn"),
               true
           )),
    _err(schema.addField(
             afw::table::Field<FluxErr>(
                 prefix + ".err", "1-sigma uncertainty on " + prefix + ".value", "dn"
             ),
             true
         ))
{}

void FluxAlgorithmMapper::apply(afw::table::BaseRecord & record, FluxAlgorithmResult const & result) {
    record.set(_value, result.value);
    record.set(_err, result.err);
    record.set(_flag, false);
}

CentroidAlgorithmMapper::CentroidAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix
) : BaseAlgorithmMapper(schema, prefix),
    _value(schema.addField(
               afw::table::Field< afw::table::Point<double> >(
                   prefix + ".value",
                   "measured centroid",
                   "pixels"
               ),
               true // replace existing fields if present
           )),
    _cov(schema.addField(
             afw::table::Field< afw::table::Covariance<afw::table::Point<float> > >(
                 prefix + ".cov",
                 "uncertainty covariance matrix for " + prefix + ".value",
                 "pixels^2"
             ),
             true // replace existing fields if present
         ))
{}

void CentroidAlgorithmMapper::apply(afw::table::BaseRecord & record, CentroidAlgorithmResult const & result) {
    record.set(_value, result.value);
    record.set(_cov, result.cov);
    record.set(_flag, false);
}

ShapeAlgorithmMapper::ShapeAlgorithmMapper(
    afw::table::Schema & schema,
    std::string const & prefix
) : BaseAlgorithmMapper(schema, prefix),
    _value(schema.addField(
               afw::table::Field< afw::table::Moments<double> >(
                   prefix + ".value",
                   "measured shape",
                   "pixels^2"
               ),
               true // replace existing fields if present
           )),
    _cov(schema.addField(
             afw::table::Field< afw::table::Covariance<afw::table::Moments<float> > >(
                 prefix + ".cov",
                 "uncertainty covariance matrix for " + prefix + ".value",
                 "pixels^4"
             ),
             true // replace existing fields if present
         ))
{}

void ShapeAlgorithmMapper::apply(afw::table::BaseRecord & record, ShapeAlgorithmResult const & result) {
    record.set(_value, result.value);
    record.set(_cov, result.cov);
    record.set(_flag, false);
}

}}} // namespace lsst::meas::base
