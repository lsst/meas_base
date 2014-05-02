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

#ifndef LSST_MEAS_BASE_Classification_h_INCLUDED
#define LSST_MEAS_BASE_Classification_h_INCLUDED

/**
 *  @file lsst/meas/base/Classification.h
 *  Implememention of the Classification algorithm under the new measurement framework
 *  Code was taken from lsst::afw::meas::algorithms with no algorithmic changes
 */

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Inputs.h"
#include "lsst/meas/base/ResultMappers.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief A C++ control class to handle ClassificationAlgorithm's configuration
 *
 *  Control fields for controlling the Classification Algorithm.
 */
class ClassificationControl {
public:
    LSST_CONTROL_FIELD(fluxRatio, double, "critical ratio of model to psf flux");
    LSST_CONTROL_FIELD(modelErrFactor, double, "correction factor for modelFlux error");
    LSST_CONTROL_FIELD(psfErrFactor, double, "correction factor for psfFlux error");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    ClassificationControl() : fluxRatio(0.925), modelErrFactor(0.0), psfErrFactor(0.0) {};
};


/// An Input struct for Classification to feed in the psfFlux and modelFlux values
//  and errors from the record object
struct ClassificationInput {

    typedef std::vector<ClassificationInput> Vector;

    float psfFlux_;
    float psfFluxSigma_;
    float modelFlux_;
    float modelFluxSigma_;

    explicit ClassificationInput(
        double psfFlux, float psfFluxSigma
    ) : psfFlux_(psfFlux), psfFluxSigma_(psfFluxSigma) {};

    ClassificationInput(afw::table::SourceRecord const & record);

    static Vector makeVector(afw::table::SourceCatalog const & catalog);

};

/**
 *  @brief Additional results for ClassificationAlgorithm
 *
 *  Unlike PsfFlux, some of Classification's outputs aren't handled by the standard FluxComponent,
 *  CentroidComponent, and ShapeComponent classes, so we have to define a Component class here
 *  to handle just those output which require special handling.
 *
 *  A corresponding ComponentMapper class is also required (see below).
 */

class ClassificationExtras {
public:
    ClassificationExtras();
    double probability;   //probability of being extended
};


/**
 *  @brief Object that transfers additional SdssShapeAlgorithm results to afw::table records
 *
 *  Because we have custom outputs, we also have to define how to transfer those outputs to
 *  records.  We just follow the pattern established by the other ComponentMapper classes.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class ClassificationExtrasMapper {
public:

    /**
     *  @brief Allocate fields in the schema and save keys for future use.
     *
     *  Unlike the standard ComponentMappers, ClassificationExtrasMappers takes a Control instance
     *  as its third argument.  It doesn't actually need it, but this is a good pattern to
     *  establish as some algorithms' outputs *will* depend on the Control object's values, and
     *  the mapper needs to have some kind of third argument in order to work with
     *  @ref measBaseResultMapperTemplates.
     *
     *  All fields should start with the given prefix and an underscore.
     */
    ClassificationExtrasMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        ClassificationControl const & control=ClassificationControl()
    );

    /// Transfer values from the result struct to the record.
    void apply(afw::table::BaseRecord & record, ClassificationExtras const & result) const;

private:
    double _fluxRatio;
    double _modelErrFactor;
    double _psfErrFactor;
    afw::table::Key<double> _probability;
};

/**
 *  @brief A measurement algorithm classifies the object by how probable it is that it is a star
 *  The ratio of the model flux and psf flux are compared to a ratio.
 */
class ClassificationAlgorithm {
public:

    /**
     *  @brief Flag bits to be used with the 'flags' data member of the Result object.
     *
     *  Inspect getFlagDefinitions() for more detailed explanations of each flag.
     *
     *  Note that we've included a final N_FLAGS value that isn't a valid flag; this is a common C++
     *  idiom for automatically counting the number of enum values, and it's required for Algorithms
     *  as the N_FLAGS value is used by the Result and ResultMapper objects.
     */
    enum FlagBits {
        N_FLAGS
    };

    /**
     *  @brief Return an array of (name, doc) tuples that describes the flags and sets the names used
     *         in catalog schemas.
     */
    static boost::array<FlagDef,N_FLAGS> const & getFlagDefinitions() {
        static boost::array<FlagDef,N_FLAGS> const flagDefs = {};
        return flagDefs;
    }

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef ClassificationControl Control;

    /**
     *  This is the type returned by apply().
     */
    typedef Result1<
        ClassificationAlgorithm,
        ClassificationExtras
    > Result;

    /// @copydoc PsfFluxAlgorithm::ResultMapper
    typedef ResultMapper1<
        ClassificationAlgorithm,
        ClassificationExtrasMapper
    > ResultMapper;


    /**
     *  In the actual overload of apply() used by the Plugin system, this is the only argument besides the
     *  Exposure being measured.  Classification needs psf and model fluxes and their errors
     */
    typedef ClassificationInput Input; // type passed to apply in addition to Exposure.

    /**
     *  @brief Create an object that transfers Result values to a record associated with the given schema
     */
    static ResultMapper makeResultMapper(
        afw::table::Schema & schema,
        std::string const & prefix,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Measure the probability of a source using the Classification algorithm.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        float psfFlux,
        float psfFluxErr,
        float modelFlux,
        float modelFluxErr,
        Control const & ctrl=Control()
    );

    /**
     *  @brief Apply the Classification to a single source using the Plugin API.
     */
    template <typename T>
    static Result apply(
        afw::image::Exposure<T> const & exposure,
        Input const & inputs,
        Control const & ctrl=Control()
    );
};

}}} // namespace lsst::meas::base

#endif // !LSST_MEAS_BASE_Classification_h_INCLUDED
