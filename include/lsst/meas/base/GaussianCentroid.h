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

#ifndef LSST_MEAS_BASE_GaussianCentroid_h_INCLUDED
#define LSST_MEAS_BASE_GaussianCentroid_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace base {

struct FittedModel {
    enum { PEAK = 0, SKY, X0, Y0, SIGMA, NPARAM };

    enum {
        BAD_GUESS = -11,
        TOO_FEW = -12,
        CHI_SQUARED = -13,
        RANGE = -14,
        BAD_WIDTH = -15,
        LOST = -16,
        DIAGONAL = -17,
        BAD_A = -18,
        CONVERGE = 1,
        ITERATE = 2,
        ALMOST = 3,
        POOR = 4
    };

    FittedModel(int status_, std::vector<double> params_, int iter_=0, double flamd_=0, double chnew_=0) :
        status(status_), params(params_), iter(iter_), flamd(flamd_), chnew(chnew_) { }
    int status;
    std::vector<double> params;
    int iter;
    double flamd;
    double chnew;
};


/**
 *  @brief A C++ control class to handle GaussianCentroidAlgorithm's configuration
 *
 *  At present, GaussianCentroidAlgorithm has no configuration options.
 *  Also, it does not currently set an error, though it should
 */
class GaussianCentroidControl {
public:
    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    GaussianCentroidControl() {}
};

/**
 *  @brief A class that calculates a centroid by fitting a circular Gaussian to the image.
 */
class GaussianCentroidAlgorithm : public SimpleAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     */
    enum {
        FAILURE=FlagHandler::FAILURE,
        NO_PEAK,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef GaussianCentroidControl Control;

    GaussianCentroidAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

private:

    // These are private so they doesn't shadow the other overloads in base classes;
    // we can still call it via the public method on the base class.  We could have
    // used a using declaration instead, but Swig had trouble with that here.

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        MeasurementError * error=NULL
    ) const;

    Control _ctrl;
    CentroidResultKey _centroidKey;
    FlagHandler _flagHandler;
    SafeCentroidExtractor _centroidExtractor;
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_GaussianCentroid_h_INCLUDED
