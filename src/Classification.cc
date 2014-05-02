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
#include "ndarray/eigen.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include "lsst/afw/table/Schema.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/base/Classification.h"


namespace lsst { namespace meas { namespace base {


ClassificationExtras::ClassificationExtras() : probability(std::numeric_limits<double>::quiet_NaN()) {}; ///< Constructor; initializes everything to NaN

ClassificationExtrasMapper::ClassificationExtrasMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ClassificationControl const & control
    ) : _fluxRatio(control.fluxRatio), _modelErrFactor(control.modelErrFactor), _psfErrFactor(control.psfErrFactor),
    _probability(
        schema.addField(
            afw::table::Field<double>(
                prefix + "_probability", "probability that the object is extended", ""
            ), true
        )
    )
{};

void ClassificationExtrasMapper::apply(afw::table::BaseRecord & record, ClassificationExtras const & result) const
{
    record.set(_probability, result.probability);
};

ClassificationAlgorithm::ResultMapper ClassificationAlgorithm::makeResultMapper(
    afw::table::Schema & schema, std::string const & name, Control const & ctrl)
{
    return ResultMapper(schema, name, ctrl);
}

template <typename T>
ClassificationAlgorithm::Result ClassificationAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    float psfFlux,
    float psfFluxErr,
    float modelFlux,
    float modelFluxErr,
    Control const & ctrl
) {
    typedef typename afw::image::Exposure<T>::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image ImageT;
    Result result;

    // set the probability by comparing the corrected model flux to the corrected psf flux
    result.probability = (ctrl.fluxRatio*modelFlux + ctrl.modelErrFactor*modelFluxErr)
        < (psfFlux + ctrl.psfErrFactor*psfFluxErr) ? 0.0 : 1.0;

    return result;
}

ClassificationInput::Vector ClassificationInput::makeVector(afw::table::SourceCatalog const & catalog) {
    Vector r;
    r.reserve(catalog.size());
    for (afw::table::SourceCatalog::const_iterator i = catalog.begin(), end=catalog.end(); i != end; ++i) {
        r.push_back(ClassificationInput(*i));
    }
    return r;
}

ClassificationInput::ClassificationInput(afw::table::SourceRecord const & record) {
    psfFlux_ = record.getPsfFlux();
    psfFluxSigma_ = record.getPsfFluxErr();
    modelFlux_ = record.getModelFlux();
    modelFluxSigma_ = record.getModelFluxErr();
}

template <typename T>
ClassificationAlgorithm::Result ClassificationAlgorithm::apply(
    afw::image::Exposure<T> const & exposure,
    Input const & inputs,
    Control const & ctrl
) {
    return apply(exposure, inputs.psfFlux_, inputs.psfFluxSigma_, inputs.modelFlux_, inputs.modelFluxSigma_, ctrl);
}

#define INSTANTIATE(T)                                                  \
    template ClassificationAlgorithm::Result ClassificationAlgorithm::apply(          \
        afw::image::Exposure<T> const & exposure,                       \
        float psfFlux,                                                  \
        float psfFluxSigma,                                             \
        float modelFlux,                                                  \
        float modelFluxSigma,                                             \
        Control const & ctrl                                            \
    );                                                                  \
    template                                                            \
    ClassificationAlgorithm::Result ClassificationAlgorithm::apply(                   \
        afw::image::Exposure<T> const & exposure,                       \
        Input const & inputs,                                           \
        Control const & ctrl                                            \
    )

INSTANTIATE(float);
INSTANTIATE(double);

}}} // namespace lsst::meas::base


