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

#include <cmath>

#include "lsst/geom/Angle.h"
#include "lsst/geom/Point.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/afw/table/BaseRecord.h"

namespace lsst {
namespace meas {
namespace base {

CentroidResult::CentroidResult()
        : x(std::numeric_limits<CentroidElement>::quiet_NaN()),
          y(std::numeric_limits<CentroidElement>::quiet_NaN()),
          xErr(std::numeric_limits<ErrElement>::quiet_NaN()),
          yErr(std::numeric_limits<ErrElement>::quiet_NaN()),
          x_y_Cov(std::numeric_limits<ErrElement>::quiet_NaN()) {}

Centroid const CentroidResult::getCentroid() const { return Centroid(x, y); }

void CentroidResult::setCentroid(Centroid const &centroid) {
    x = centroid.getX();
    y = centroid.getY();
}

CentroidCov const CentroidResult::getCentroidErr() const {
    CentroidCov m;
    m << xErr * xErr, x_y_Cov, x_y_Cov, yErr * yErr;
    return m;
}

void CentroidResult::setCentroidErr(CentroidCov const &matrix) {
    xErr = std::sqrt(matrix(0, 0));
    yErr = std::sqrt(matrix(1, 1));
    x_y_Cov = matrix(0, 1);
}

void CentroidResult::setCentroidErr(ErrElement _xErr, ErrElement _yErr) {
    xErr = _xErr;
    yErr = _yErr;
    x_y_Cov = 0.0;
}

CentroidResultKey CentroidResultKey::addFields(afw::table::Schema &schema, std::string const &name,
                                               std::string const &doc, UncertaintyEnum uncertainty) {
    CentroidResultKey r;
    r._centroid = afw::table::PointKey<CentroidElement>::addFields(schema, name, doc, "pixel");
    if (uncertainty != NO_UNCERTAINTY) {
        std::vector<afw::table::Key<ErrElement> > sigma(2);
        std::vector<afw::table::Key<ErrElement> > cov;
        sigma[0] = schema.addField<ErrElement>(schema.join(name, "xErr"),
                                               "1-sigma uncertainty on x position", "pixel");
        sigma[1] = schema.addField<ErrElement>(schema.join(name, "yErr"),
                                               "1-sigma uncertainty on y position", "pixel");
        if (uncertainty == FULL_COVARIANCE) {
            cov.push_back(schema.addField<ErrElement>(schema.join(name, "x_y_Cov"),
                                                      "uncertainty covariance in x and y", "pixel^2"));
        }
        r._centroidErr = afw::table::CovarianceMatrixKey<ErrElement, 2>(sigma, cov);
    }
    return r;
}

namespace {

std::vector<std::string> getNameVector() {
    std::vector<std::string> v;
    v.push_back("x");
    v.push_back("y");
    return v;
}

}  // namespace

CentroidResultKey::CentroidResultKey(afw::table::SubSchema const &s) : _centroid(s) {
    static std::vector<std::string> names = getNameVector();  // C++11 TODO: just use initializer list
    try {
        _centroidErr = afw::table::CovarianceMatrixKey<ErrElement, 2>(s, names);
    } catch (pex::exceptions::NotFoundError &) {
    }
}

CentroidResult CentroidResultKey::get(afw::table::BaseRecord const &record) const {
    CentroidResult r;
    r.setCentroid(record.get(_centroid));
    if (_centroidErr.isValid()) {
        r.setCentroidErr(record.get(_centroidErr));
    }
    return r;
}

void CentroidResultKey::set(afw::table::BaseRecord &record, CentroidResult const &value) const {
    record.set(_centroid, value.getCentroid());
    if (_centroidErr.isValid()) {
        record.set(_centroidErr, value.getCentroidErr());
    }
}

CentroidTransform::CentroidTransform(std::string const &name, afw::table::SchemaMapper &mapper)
        : BaseTransform{name} {
    // Map the flag through to the output
    mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(name + "_flag").key);

    // Add keys for the coordinates
    auto &s = mapper.editOutputSchema();
    _coordKey = afw::table::CoordKey::addFields(s, name, "ICRS coordinates");

    // If the centroid has an error key we also include one on the celestial
    // coordinates; otherwise, it isn't necessary. Note that if we provide for
    // errors in celestial coordinates, we always need the full covariance.
    if (CentroidResultKey(mapper.getInputSchema()[name]).getCentroidErr().isValid()) {
        std::vector<afw::table::Key<ErrElement> > sigma(2);
        std::vector<afw::table::Key<ErrElement> > cov(1);
        sigma[0] = s.addField<ErrElement>(s.join(name, "raErr"), "1-sigma uncertainty on RA", "rad");
        sigma[1] = s.addField<ErrElement>(s.join(name, "decErr"), "1-sigma uncertainty on dec", "rad");
        cov[0] = s.addField<ErrElement>(s.join(name, "ra_dec_Cov"), "Uncertainty covariance in RA and dec",
                                        "rad^2");
        _coordErrKey = afw::table::CovarianceMatrixKey<ErrElement, 2>(sigma, cov);
    }
}

void CentroidTransform::operator()(afw::table::SourceCatalog const &inputCatalog,
                                   afw::table::BaseCatalog &outputCatalog, afw::geom::SkyWcs const &wcs,
                                   afw::image::Calib const &calib) const {
    checkCatalogSize(inputCatalog, outputCatalog);
    CentroidResultKey centroidResultKey(inputCatalog.getSchema()[_name]);

    afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
    afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();

    for (; inSrc != inputCatalog.end() && outSrc != outputCatalog.end(); ++inSrc, ++outSrc) {
        CentroidResult centroidResult = centroidResultKey.get(*inSrc);

        _coordKey.set(*outSrc, wcs.pixelToSky(centroidResult.getCentroid()));

        if (centroidResultKey.getCentroidErr().isValid()) {
            CentroidCov centroidCov = centroidResult.getCentroidErr();
            if (!(std::isnan(centroidCov(0, 0)) || std::isnan(centroidCov(1, 1)))) {
                auto transform = wcs.linearizePixelToSky(centroidResult.getCentroid(), geom::radians)
                                         .getLinear()
                                         .getMatrix();
                _coordErrKey.set(*outSrc, (transform * centroidResult.getCentroidErr().cast<double>() *
                                           transform.transpose())
                                                  .cast<ErrElement>());
            }
        }
    }
}

// Add a key "flag_resetToPeak" if something is wrong with the centroid
// Save a key to this algorithm's Centroid, as well as the general failure flag
CentroidChecker::CentroidChecker(afw::table::Schema &schema, std::string const &name, bool doFootprintCheck,
                                 double maxDistFromPeak)
        : _doFootprintCheck(doFootprintCheck), _maxDistFromPeak(maxDistFromPeak) {
    _resetKey = schema.addField<afw::table::Flag>(schema.join(name, "flag_resetToPeak"),
                                                  "set if CentroidChecker reset the centroid");
    _failureKey = schema.find<afw::table::Flag>(schema.join(name, "flag")).key;
    _xKey = schema.find<CentroidElement>(schema.join(name, "x")).key;
    _yKey = schema.find<CentroidElement>(schema.join(name, "y")).key;
}

//  Set the centroid to the first footprint if the centroid is eithe more than _maxDistFromPeak
//  pixels from the centroid, or if it is outside the footprint.
bool CentroidChecker::operator()(afw::table::SourceRecord &record) const {
    CentroidElement x = record.get(_xKey);
    CentroidElement y = record.get(_yKey);

    if (!_doFootprintCheck && _maxDistFromPeak < 0.0) {
        return false;
    }
    PTR(afw::detection::Footprint) footprint = record.getFootprint();
    if (!footprint) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, "No Footprint attached to record");
    }
    if (footprint->getPeaks().empty()) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, "Footprint has no peaks; cannot verify centroid.");
    }
    CentroidElement footX = footprint->getPeaks().front().getFx();
    CentroidElement footY = footprint->getPeaks().front().getFy();
    double distsq = (x - footX) * (x - footX) + (y - footY) * (y - footY);
    if ((_doFootprintCheck && !footprint->contains(geom::Point2I(geom::Point2D(x, y)))) ||
        ((_maxDistFromPeak > 0) && (distsq > _maxDistFromPeak * _maxDistFromPeak))) {
        record.set(_xKey, footX);
        record.set(_yKey, footY);
        record.set(_failureKey, true);
        record.set(_resetKey, true);
        return true;
    }
    return false;
}
}  // namespace base
}  // namespace meas
}  // namespace lsst
