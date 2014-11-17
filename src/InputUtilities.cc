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

#include "lsst/utils/ieee.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace base {

SafeCentroidExtractor::SafeCentroidExtractor(
    afw::table::Schema & schema,
    std::string const & name,
    bool isCentroider
) :
    _name(name),
    _isCentroider(isCentroider)
{
    // Instead of aliasing e.g. MyAlgorithm_flag_badCentroid->slot_Centroid_flag, we actually
    // look up the target of slot_Centroid_flag, and alias that to MyAlgorithm_flag_badCentroid.
    // That way, if someone changes the slots later, after we've already done the measurement,
    // this alias still points to the right thing.
    std::string slotFlagName = schema.getAliasMap()->apply(schema.join("slot", "Centroid", "flag"));
    if (_isCentroider) {
        if (slotFlagName != schema.join(name, "flag")) {
            // only setup the alias if this isn't the slot algorithm itself (otherwise it'd be circular)
            schema.getAliasMap()->set(schema.join(name, "flag", "badInitialCentroid"), slotFlagName);
        }
    } else {
        schema.getAliasMap()->set(schema.join(name, "flag", "badCentroid"), slotFlagName);
    }
}

namespace {

afw::geom::Point2D extractPeak(afw::table::SourceRecord const & record, std::string const & name) {
    afw::geom::Point2D result;
    PTR(afw::detection::Footprint) footprint = record.getFootprint();
    if (!footprint) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeError,
            (boost::format("%s: Centroid slot value is NaN, but no Footprint attached to record")
             % name).str()
        );
    }
    if (footprint->getPeaks().empty()) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeError,
            (boost::format("%s: Centroid slot value is NaN, but Footprint has no Peaks")
             % name).str()
        );
    }
    result.setX(footprint->getPeaks().front()->getFx());
    result.setY(footprint->getPeaks().front()->getFy());
    return result;
}

} // anonymous

afw::geom::Point2D SafeCentroidExtractor::operator()(
    afw::table::SourceRecord & record,
    FlagHandler const & flags
) const {
    if (!record.getTable()->getCentroidKey().isValid()) {
        if (_isCentroider) {
            return extractPeak(record, _name);
        } else {
            throw LSST_EXCEPT(
                FatalAlgorithmError,
                (boost::format("%s requires a centroid, but the centroid slot is not defined") % _name).str()
            );
        }
    }
    afw::geom::Point2D result = record.getCentroid();
    if (utils::isnan(result.getX()) || utils::isnan(result.getY())) {
        if (!record.getTable()->getCentroidFlagKey().isValid()) {
            if (_isCentroider) {
                return extractPeak(record, _name);
            } else {
                throw LSST_EXCEPT(
                    pex::exceptions::RuntimeError,
                    (boost::format("%s: Centroid slot value is NaN, but there is no Centroid slot flag "
                                   "(is the executionOrder for %s lower than that of the slot Centroid?)")
                     % _name % _name).str()
                );
            }
        }
        if (!record.getCentroidFlag() && !_isCentroider) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                (boost::format("%s: Centroid slot value is NaN, but the Centroid slot flag is not set "
                               "(is the executionOrder for %s lower than that of the slot Centroid?)")
                 % _name % _name).str()
            );
        }
        result = extractPeak(record, _name);
        if (!_isCentroider) {
            // set the general flag, because using the Peak might affect the current measurement
            flags.setValue(record, FlagHandler::FAILURE, true);
        }
    } else if (!_isCentroider && record.getTable()->getCentroidFlagKey().isValid()
               && record.getCentroidFlag()) {
        // we got a usable value, but the centroid flag is still be set, and that might affect
        // the current measurement
        flags.setValue(record, FlagHandler::FAILURE, true);
    }
    return result;
}

SafeShapeExtractor::SafeShapeExtractor(afw::table::Schema & schema, std::string const & name) :
    _name(name)
{
    // Instead of aliasing e.g. MyAlgorithm_flag_badShape->slot_Shape_flag, we actually
    // look up the target of slot_Shape_flag, and alias that to MyAlgorithm_flag_badCentroid.
    // That way, if someone changes the slots later, after we've already done the measurement,
    // this alias still points to the right thing.
    std::string slotFlagName = schema.getAliasMap()->apply(schema.join("slot", "Shape", "flag"));
    schema.getAliasMap()->set(schema.join(name, "flag", "badShape"), slotFlagName);
}

afw::geom::ellipses::Quadrupole SafeShapeExtractor::operator()(
    afw::table::SourceRecord & record,
    FlagHandler const & flags
) const {
    if (!record.getTable()->getShapeKey().isValid()) {
        throw LSST_EXCEPT(
            FatalAlgorithmError,
            (boost::format("%s requires a shape, but the shape slot is not defined") % _name).str()
        );
    }
    afw::geom::ellipses::Quadrupole result = record.getShape();
    if (utils::isnan(result.getIxx()) || utils::isnan(result.getIyy()) || utils::isnan(result.getIxy())) {
        if (!record.getTable()->getShapeFlagKey().isValid()) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                (boost::format("%s: Shape slot value is NaN, but there is no Shape slot flag "
                               "(is the executionOrder for %s lower than that of the slot Shape?)")
                 % _name % _name).str()
            );
        }
        if (!record.getShapeFlag()) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                (boost::format("%s: Shape slot value is NaN, but the Shape slot flag is not set "
                               "(is the executionOrder for %s lower than that of the slot Shape?)")
                 % _name % _name).str()
            );
        }
        throw LSST_EXCEPT(
            MeasurementError,
            (boost::format("%s: Shape needed, and Shape slot measurement failed.") % _name).str(),
            FlagHandler::FAILURE
        );
    } else if (record.getTable()->getShapeFlagKey().isValid() && record.getShapeFlag()) {
        // we got a usable value, but the shape flag might still be set, and that might affect
        // the current measurement
        flags.setValue(record, FlagHandler::FAILURE, true);
    }
    return result;
}

}}} // lsst::meas::base
