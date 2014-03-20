// -*- LSST-C++ -*-
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

#include "lsst/utils/ieee.h"
#include "lsst/meas/algorithms/CorrectFluxes.h"
#include "lsst/afw/detection/Psf.h"

namespace lsst { namespace meas { namespace algorithms {


typedef std::pair<afw::table::KeyTuple<afw::table::Flux>,ScaledFlux::KeyTuple> ScaledFluxKeys;

// Helper function to get a value from a record and return a default value when the key is invalid.
template <typename T>
inline typename afw::table::Field<T>::Value get(
    afw::table::SourceRecord const & record,
    afw::table::Key<T> const & key,
    typename afw::table::Field<T>::Value default_ = false
) {
    return key.isValid() ? record.get(key) : default_;
}

// Helper functor to set the general flag for a flux measurement
struct SetMainFlag {

    explicit SetMainFlag(afw::table::SourceRecord * record_) : record(record_) {}

    void operator()(ScaledFluxKeys const & keys) const {
        if (keys.first.flag.isValid()) {
            record->set(keys.first.flag, true);
        }
    }
    
    afw::table::SourceRecord * record;
};

// Helper functor to apply the aperture correction to a flux measurement
struct ApplyApCorr {

    ApplyApCorr(double apCorr_, afw::table::SourceRecord * record_) : apCorr(apCorr_), record(record_) {}

    void operator()(ScaledFluxKeys const & keys) const {
        (*record)[keys.first.meas] *= apCorr;
        if (keys.first.err.isValid()) {
            (*record)[keys.first.err] *= apCorr;
        }
    }
    
    double apCorr;
    afw::table::SourceRecord * record;
};


class CorrectFluxes : public Algorithm {
public:

    CorrectFluxes(
        CorrectFluxesControl const & ctrl,
        afw::table::Schema & schema,
        AlgorithmMap const & others
    );

private:

    void _apply(
        afw::table::SourceRecord & source,
        PTR(afw::detection::Psf const) psf,
        afw::geom::Point2D const & center
    ) const;
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const {
        return this->_apply(source, exposure.getPsf(), center);
    }

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(CorrectFluxes);

    afw::table::Key<float> _apCorrKey;
    afw::table::Key<afw::table::Flag> _apCorrFlagKey;
    ScaledFluxKeys _canonical;
    std::vector<ScaledFluxKeys> _others;
};

CorrectFluxes::CorrectFluxes(
    CorrectFluxesControl const & ctrl,
    afw::table::Schema & schema,
    AlgorithmMap const & others
) : Algorithm(ctrl) {

    if (ctrl.doApCorr) {
        _apCorrKey = schema.addField<float>(
            ctrl.name + ".apcorr",
            (boost::format("correction applied to model fluxes to match to %f-pixel aperture")
             % ctrl.apCorrRadius).str()
        );
        _apCorrFlagKey = schema.addField<afw::table::Flag>(
            ctrl.name + ".apcorr.flags",
            "flag set if aperture correction failed"
        );
    }

    for (AlgorithmMap::const_iterator i = others.begin(); i != others.end(); ++i) {
        CONST_PTR(ScaledFlux) asScaledFlux = boost::dynamic_pointer_cast<ScaledFlux const>(i->second);
        if (i->second->getControl().name == ctrl.canonicalFluxName) {
            if (!asScaledFlux) {
                throw LSST_EXCEPT(
                    pex::exceptions::InvalidParameterException,
                    (boost::format("Canonical flux (%s) is not an instance of ScaledFlux")
                     % ctrl.canonicalFluxName).str()
                );
            }
            int fluxCount = asScaledFlux->getFluxCount();
            if (ctrl.canonicalFluxIndex < 0 || ctrl.canonicalFluxIndex >= fluxCount) {
                throw LSST_EXCEPT(
                    pex::exceptions::InvalidParameterException,
                    (boost::format("Invalid index (%d) for canonical flux (must be between 0 and %d)")
                     % ctrl.canonicalFluxIndex % (fluxCount - 1)).str()
                );
            }
            for (int i = 0; i < fluxCount; ++i) {
                if (i == ctrl.canonicalFluxIndex) {
                    _canonical =
                        ScaledFluxKeys(
                            asScaledFlux->getFluxKeys(i),
                            asScaledFlux->getFluxCorrectionKeys(i)
                        );
                } else {
                    _others.push_back(
                        ScaledFluxKeys(
                            asScaledFlux->getFluxKeys(i),
                            asScaledFlux->getFluxCorrectionKeys(i)
                        )
                    );
                }
            }
        } else if (asScaledFlux) {
            int fluxCount = asScaledFlux->getFluxCount();
            for (int i = 0; i < fluxCount; ++i) {
                _others.push_back(
                    ScaledFluxKeys(
                        asScaledFlux->getFluxKeys(i),
                        asScaledFlux->getFluxCorrectionKeys(i)
                    )
                );
            }
        }
    } // for i in others

    if (ctrl.doTieScaledFluxes && !_canonical.first.meas.isValid()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot tie scaled fluxes without a canonical flux measurement."
        );
    }
}



}}} // namespace lsst::meas::algorithms
