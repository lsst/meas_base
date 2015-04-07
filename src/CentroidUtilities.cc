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

#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/afw/table/BaseRecord.h"

namespace lsst { namespace meas { namespace base {

CentroidResult::CentroidResult() :
    x(std::numeric_limits<CentroidElement>::quiet_NaN()),
    y(std::numeric_limits<CentroidElement>::quiet_NaN()),
    xSigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    ySigma(std::numeric_limits<ErrElement>::quiet_NaN()),
    x_y_Cov(std::numeric_limits<ErrElement>::quiet_NaN())
{}

Centroid const CentroidResult::getCentroid() const { return Centroid(x, y); }

void CentroidResult::setCentroid(Centroid const & centroid) {
    x = centroid.getX();
    y = centroid.getY();
}

CentroidCov const CentroidResult::getCentroidErr() const {
    CentroidCov m;
    m <<
        xSigma*xSigma, x_y_Cov,
        x_y_Cov, ySigma*ySigma;
    return m;
}

void CentroidResult::setCentroidErr(CentroidCov const & matrix) {
    xSigma = std::sqrt(matrix(0, 0));
    ySigma = std::sqrt(matrix(1, 1));
    x_y_Cov = matrix(0, 1);
}

void CentroidResult::setCentroidErr(ErrElement _xSigma, ErrElement _ySigma) {
    xSigma = _xSigma;
    ySigma = _ySigma;
    x_y_Cov = 0.0;
}

CentroidResultKey CentroidResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name,
    std::string const & doc,
    UncertaintyEnum uncertainty
) {
    CentroidResultKey r;
    r._centroid = afw::table::PointKey<CentroidElement>::addFields(
        schema,
        name,
        doc,
        "pixels"
    );
    if (uncertainty != NO_UNCERTAINTY) {
        std::vector< afw::table::Key<ErrElement> > sigma(2);
        std::vector< afw::table::Key<ErrElement> > cov;
        sigma[0] = schema.addField<ErrElement>(
            schema.join(name, "xSigma"), "1-sigma uncertainty on x position", "pixels"
        );
        sigma[1] = schema.addField<ErrElement>(
            schema.join(name, "ySigma"), "1-sigma uncertainty on y position", "pixels"
        );
        if (uncertainty == FULL_COVARIANCE) {
            cov.push_back(
                schema.addField<ErrElement>(
                    schema.join(name, "x_y_Cov"), "uncertainty covariance in x and y", "pixels^2"
                )
            );
        }
        r._centroidErr = afw::table::CovarianceMatrixKey<ErrElement,2>(sigma, cov);
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

} // anonymous

CentroidResultKey::CentroidResultKey(afw::table::SubSchema const & s) :
    _centroid(s)
{
    static std::vector<std::string> names = getNameVector(); // C++11 TODO: just use initializer list
    try {
        _centroidErr = afw::table::CovarianceMatrixKey<ErrElement,2>(s, names);
    } catch (pex::exceptions::NotFoundError &) {}
}

CentroidResult CentroidResultKey::get(afw::table::BaseRecord const & record) const {
    CentroidResult r;
    r.setCentroid(record.get(_centroid));
    if (_centroidErr.isValid()) {
        r.setCentroidErr(record.get(_centroidErr));
    }
    return r;
}

void CentroidResultKey::set(afw::table::BaseRecord & record, CentroidResult const & value) const {
    record.set(_centroid, value.getCentroid());
    if (_centroidErr.isValid()) {
        record.set(_centroidErr, value.getCentroidErr());
    }
}

CentroidTransform::CentroidTransform(
    std::string const & name,
    afw::table::SchemaMapper & mapper
) :
    BaseTransform{name}
{
    // Map the flag through to the output
    mapper.addMapping(mapper.getInputSchema().find<afw::table::Flag>(name + "_flag").key);

    // Add keys for the coordinates
    auto & s = mapper.editOutputSchema();
    _coordKey = afw::table::CoordKey::addFields(s, name, "ICRS coordinates");

    // If the centroid has an error key we also include one on the celestial
    // coordinates; otherwise, it isn't necessary. Note that if we provide for
    // errors in celestial coordinates, we always need the full covariance.
    if (CentroidResultKey(mapper.getInputSchema()[name]).getCentroidErr().isValid()) {
        std::vector< afw::table::Key<ErrElement> > sigma(2);
        std::vector< afw::table::Key<ErrElement> > cov(1);
        sigma[0] = s.addField<ErrElement>(s.join(name, "raSigma"), "Uncertainty on RA", "radians");
        sigma[1] = s.addField<ErrElement>(s.join(name, "decSigma"), "Uncertainty on dec", "radians");
        cov[0] = s.addField<ErrElement>(s.join(name, "ra_dec_Cov"),
                                             "Uncertainty covariance in RA and dec", "radians^2");
        _coordErrKey = afw::table::CovarianceMatrixKey<ErrElement,2>(sigma, cov);
    }

}

void CentroidTransform::operator()(
    afw::table::SourceCatalog const & inputCatalog,
    afw::table::BaseCatalog & outputCatalog,
    afw::image::Wcs const & wcs,
    afw::image::Calib const & calib
) const {
    checkCatalogSize(inputCatalog, outputCatalog);
    CentroidResultKey centroidResultKey(inputCatalog.getSchema()[_name]);

    afw::table::SourceCatalog::const_iterator inSrc = inputCatalog.begin();
    afw::table::BaseCatalog::iterator outSrc = outputCatalog.begin();

    for (; inSrc < inputCatalog.end() && outSrc < outputCatalog.end(); ++inSrc, ++outSrc) {
        CentroidResult centroidResult = centroidResultKey.get(*inSrc);

        _coordKey.set(*outSrc, *wcs.pixelToSky(centroidResult.getCentroid()));

        if (centroidResultKey.getCentroidErr().isValid()) {
            CentroidCov centroidCov = centroidResult.getCentroidErr();
            if (!(utils::isnan(centroidCov(0,0)) || utils::isnan(centroidCov(1,1)))) {
                auto transform = wcs.linearizePixelToSky(centroidResult.getCentroid(),
                                                         afw::geom::radians).getLinear().getMatrix();
                _coordErrKey.set(*outSrc, (transform * centroidResult.getCentroidErr().cast<double>() *
                                           transform.transpose()).cast<ErrElement>());
            }
        }
    }
}

}}} // lsst::meas::base
