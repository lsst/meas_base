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

#ifndef LSST_MEAS_BASE_Results_h_INCLUDED
#define LSST_MEAS_BASE_Results_h_INCLUDED

#include <bitset>

#include "boost/array.hpp"
#include "Eigen/Core"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"

namespace lsst { namespace meas { namespace base {

/**
 *  @brief Exception to be thrown when a measurement algorithm experiences a known failure mode.
 *
 *  In addition to the usual message, MeasurementError must be constructed with the bit of the
 *  algorithm-specific flag that indicates the known failure mode.  This allows the measurement
 *  framework to set that flag upon failure.  Typically, this flag bit is also used to look up
 *  the message from the algorithm classes FlagDefinition list; the common pattern is:
 *  @code
 *      if (badThingHappened) {
 *          // BAD_THING is an enum value from the Algorithm's FlagBits enum; getFlagDefinitions()
 *          // is a static method algorithms are expected to define.
 *          throw LSST_EXCEPT(MeasurementError, getFlagDefinitions()[BAD_THING), BAD_THING);
 *       }
 *  @endcode
 */
class MeasurementError : public pex::exceptions::RuntimeErrorException {
public:

    /// Constructor; should only be invoked by the LSST_EXCEPT macro (see class docs)
    MeasurementError(LSST_EARGS_TYPED, std::size_t flagBit) :
        pex::exceptions::RuntimeErrorException(LSST_EARGS_UNTYPED),
        _flagBit(flagBit)
    {}

    /// Return the flag bit associated with the error.
    std::size_t getFlagBit() const { return _flagBit; }

    virtual char const* getType(void) const throw() { return "lsst::meas::base::MeasurementError *"; };

    virtual lsst::pex::exceptions::Exception* clone(void) const {
        return new MeasurementError(*this);
    };

private:
    std::size_t _flagBit;
};

/**
 *  @brief Exception to be thrown when a measurement algorithm encounters a NaN or infinite pixel.
 *
 *  When caught by the plugin framework, this exception will generate a log message.
 */
LSST_EXCEPTION_TYPE(PixelValueError, lsst::pex::exceptions::DomainErrorException,
                    lsst::meas::base::PixelValueError);

/**
 *  @brief An enum used to specify how much uncertainty information measurement algorithms provide.
 *
 *  Currently, only ResultMappers (not Results) make use of these distinctions; Result structs always
 *  have data members that could hold the full-covariance, but may set some of these to NaN.  The
 *  ResultMapper knows that their individual compoents may only produce a sigma, or no error field
 *  at all.  When the Results are written to the output record, this information is heeded.
 */
enum UncertaintyEnum {
    NO_UNCERTAINTY = 0, ///< Algorithm provides no uncertainy information at all
    SIGMA_ONLY = 1,     ///< Only the diagonal elements of the covariance matrix are provided
    FULL_COVARIANCE = 2 ///< The full covariance matrix is provided
};

//@{ Typedefs that define the C++ types we typically use for common measurements
typedef float Flux;
typedef float ErrElement;
typedef double CentroidElement;
typedef double ShapeElement;
typedef afw::geom::Point<CentroidElement,2> Centroid;
typedef Eigen::Matrix<ErrElement,2,2,Eigen::DontAlign> CentroidCov;
typedef afw::geom::ellipses::Quadrupole Shape;
typedef Eigen::Matrix<ErrElement,3,3,Eigen::DontAlign> ShapeCov;
//@}

/**
 *  @brief Simple POD struct used to define and document flags
 *
 *  Each Algorithm class should define a static getFlagDefinitions() method that returns a const reference
 *  to a boost::array of FlagDef objects, and that's generally the only time this class will be used.  This
 *  array is then used by the ResultMapper classes to add flag fields to afw::table::Schema objects, and by
 *  algorithm code to retrieve error messages when throwing MeasurementError.
 *
 *  When we switch to C++11, we can make the attributes std::strings, but at least for now we'll use
    C-strings so we can create these arrays using initializer lists even in C++98.
 */
struct FlagDef {
    char const * name;
    char const * doc;
};

/**
 *  @defgroup measBaseResults Results and ResultMappers
 *
 *  In order to be wrapped as Plugins, Algorithm classes must have apply() (and optionally applyN()) methods
 *  that return a Result containing the outputs of the algorithm. A Result struct is an aggregation of
 *  one or more "Components". A Component is a single measurement, such as a position and position error,
 *  or a flux and flux error.  
 *
 *  A standard set of Components (Flux, Centroid, Shape) are defined in the lsst::meas::base namespace.
 *  A FlagsComponent is also defined for returning error flags.
 *
 *  In most cases, Algorithms will create their Result from Components by instantiating one of the 
 *  @ref measBaseResultTemplates, which can be used to aggregate one or more standard Components
 *  (FluxComponent, CentroidComponent, ShapeComponent) with a set flags (FlagsComponent).  Classes whose
 *  outputs cannot be categorized this way should define their own custom "Extras" Component.
 *  (See SdssShapeAlgorithm as an example, or @ref measBaseResultTemplates for more documentation.)
 *
 *  A ResultMapper transfers the Result values to a afw::table::SourceRecord in the Plugin wrapper.
 *  An algorithm must therefore also provide a ResultMapper typedef corresponding to its Result typedef.
 *  And it must provide a makeResultMapper static method to create its ResultMapper.
 *
 *  When the Result is one of the standard ResultN templates, the ResultMapper should be one of the standard
 *  ResultMapperN templates with the same Components.  Note that if the Result is made up of standard
 *  Components, a corresponding ComponentMapper (FluxComponentMapper, CentroidComponentMapper,
 *  ShapeComponentMapper, and FlagsComponentMapper) is supplied, and creating a ResultMapper is trivial.
 *  If the Component is not one of the Standard Componets, you will also need to provide a mapper for
 *  the "extras" (again, see SdssShapeAlgorithm) for an example.
 *
 *  @{
 */

/**
 *  @brief A reusable result struct component for flags.
 *
 *  All algorithms should include a FlagsComponent in their Result struct (all of the ResultN templates do)
 *  to provide detailed information about different failure modes.  In general, however, an Algorithm should
 *  only set flags directly in the Result struct for non-fatal errors, even if the Result struct should
 *  have bits for both fatal and non-fatal errors.  For fatal errors, the MeasurementError exception
 *  should be thrown, and the appropriate bit set should be set by the exception handler.
 *
 *  In addition, while the afw::table::Schema for each Plugin will contain a Flag that indicates a general
 *  failure of that plugin, such a Flag should NOT be defined explicitly by the Algorithm, and will not be
 *  included in the FlagsComponent; instead, the Algorithm should simply throw an exception, and the Plugin
 *  framework will ensure that the general failure bit is set in this case (in addition to any other bit
 *  attached to the MeasurementError, if that is what is thrown).
 *
 *  In order to make use of FlagsComponent, an Algorithm class must:
 *   - Define an enum FlagBits that enumerates the flag values, and includes a terminating N_FLAGS enum
 *     value that declares the total number of flags.  Note that the enum values should be the bit index,
 *     NOT an integer bitmask.
 *   - Define a static makeFlagDefinitions() method that returns a const reference to
 *     @c boost::array<FlagDef,N_FLAGS>. This contains schema-appropriate names and descriptions for
 *     each flag (to create the actual schema field name, these names will be appended to a
 *     "<PluginName>_flag_" prefix).  Schema-appropriate flag names begin with a lower-case letter
 *     and are camelCase with no spaces.
 */
template <typename Algorithm>
struct FlagsComponent {

    /// Return the flag value associated with the given bit
    bool getFlag(typename Algorithm::FlagBits bit) const { return _flags[bit]; }

    /// Set the flag value associated with the given bit
    void setFlag(typename Algorithm::FlagBits bit, bool value=true) { _flags[bit] = value; }

    /// Clear (i.e. set to false) the flag associated with the given bit
    void unsetFlag(typename Algorithm::FlagBits bit) { _flags[bit] = false; }

private:

    template <typename A> friend class FlagsComponentMapper;

    std::bitset<Algorithm::N_FLAGS> _flags;
};

/**
 *  @brief A reusable component for result structs for flux measurements.
 *
 *  Flux measurements and their errors should always be in DN.
 *
 *  @todo figure out how to handle arrays of fluxes from e.g. aperture photometry
 */
struct FluxComponent {
    Flux flux; ///< Measured flux in DN.
    ErrElement fluxSigma; ///< 1-Sigma error (sqrt of variance) on flux in DN.

    FluxComponent(); ///< Constructor; initializes everything to NaN.
};

/**
 *  @brief A reusable component for result structs for centroid or other position measurements.
 *
 *  Centroid measurements and their errors should always be in pixels, relative to the image's xy0.
 */
struct CentroidComponent {
    CentroidElement x; ///< x (column) coordinate of the measured position
    CentroidElement y; ///< y (row) coordinate of the measured position
    ErrElement xSigma; ///< 1-Sigma uncertainty on x (sqrt of variance)
    ErrElement ySigma; ///< 1-Sigma uncertainty on y (sqrt of variance)
    ErrElement x_y_Cov; ///< x,y term in the uncertainty convariance matrix

    /// Return a Point object containing the measured x and y
    Centroid const getCentroid() const { return Centroid(x, y); }

    /// Return the 2x2 symmetric covariance matrix, with rows and columns ordered (x, y)
    CentroidCov const getCov() const {
        CentroidCov m;
        m <<
            xSigma*xSigma, x_y_Cov,
            x_y_Cov, ySigma*ySigma;
        return m;
    }

    CentroidComponent(); ///< Constructor; initializes everything to NaN.

};

/**
 *  @brief A reusable component for result structs for moments-based shape measurements.
 *
 *  Shape measurements and their errors should always be in pixels coordinates.  This struct should generally
 *  be preferred over a custom struct with other ellipse parametrizations unless the measurement takes place
 *  in another parametrization and a transformation to this one would result in a loss of information or
 *  obfuscate the results of the measurement (i.e. use this one unless you have a good reason not to).
 */
struct ShapeComponent {
    ShapeElement xx; // image or model second moment for x^2
    ShapeElement yy; // image or model second moment for y^2
    ShapeElement xy; // image or model second moment for xy^2
    ErrElement xxSigma; ///< 1-Sigma uncertainty on xx (sqrt of variance)
    ErrElement yySigma; ///< 1-Sigma uncertainty on yy (sqrt of variance)
    ErrElement xySigma; ///< 1-Sigma uncertainty on xy (sqrt of variance)
    ErrElement xx_yy_Cov; ///< xx,yy term in the uncertainty convariance matrix
    ErrElement xx_xy_Cov; ///< xx,xy term in the uncertainty convariance matrix
    ErrElement yy_xy_Cov; ///< yy,xy term in the uncertainty convariance matrix

    /**
     *  @brief Return an afw::geom::ellipses object corresponding to xx, yy, xy.
     *
     *  This method can be used to return an average radius for the measured shape, e.g.
     *  @c getShape().getDeterminantRadius()
     */
    Shape const getShape() const { return Shape(xx, yy, xy); }

    /// Return the 3x3 symmetric covariance matrix, with rows and columns ordered (xx, yy, xy)
    ShapeCov const getCov() const {
        ShapeCov m;
        m <<
            xxSigma*xxSigma, xx_yy_Cov, xx_xy_Cov,
            xx_yy_Cov, yySigma*yySigma, yy_xy_Cov,
            xx_xy_Cov, yy_xy_Cov, xySigma*xySigma;
        return m;
    }

    ShapeComponent(); ///< Constructor; initializes everything to NaN.

};

/**
 *  @defgroup measBaseResultTemplates the ResultN templates
 *
 *  All Algorithm classes should declare a typedef 'Result'.  This typedef tells the Plugin wrapper system
 *  the type which will be returned by their apply() method.  A Result is typically an instantiation of
 *  one of ResultN templates, with N Components (the template arguments).  There are a few requirements
 *  and recommendations for classes that can be used as Components:
 *   - They must be default constructable, copyable, and assignable.
 *   - The default constructor should initialize all data members to sensible values (usually NaN for
 *     floating point values).
 *   - Public data members are allowed, but should be limited to scalar POD types supported by afw::table.
 *   - When possible the typedefs used in the standard Component structs should used (i.e. ErrElement for
 *     uncertainty values).
 *
 *  Custom Components must be accompanied by an associated ComponentMapper (see
 *  @ref measBaseResultMapperTemplates).  The ComponentMapper is somewhat more complicated than the
 *  Component itself.
 *
 *  Note that these templates all implicitly include a FlagsComponent.
 *
 *  @internal
 *  We use multiple inheritance here (following a pattern commonly used in implementations of std::tuple)
 *  rather than composition to keep user code concise and boilerplate-free,
 *  and also to allow the Python/C++ code for accessing a Result struct to more closely resemble the
 *  afw::table fields the Result struct data members correspond to.  Because all of the component classes
 *  are strictly nonpolymorphic, it's better to just think of this as a language-specific technique for
 *  struct concatenation that true OO multiple-inheritance.
 *
 *  In C++11, we could consider making these variadic, but the fact that we'll likely never need more than
 *  four Components means the implementation of a variadic version (let alone its Swig wrappers) would
 *  almost certainly require more code (and be harder to maintain) than the explicit version here.
 *  @endinternal
 * @{
 */
template <typename Algorithm, typename T1>
struct Result1 : public T1, public FlagsComponent<Algorithm> {};

template <typename Algorithm, typename T1, typename T2>
struct Result2 : public T1, public T2, public FlagsComponent<Algorithm> {};

template <typename Algorithm, typename T1, typename T2, typename T3>
struct Result3 : public T1, public T2, public T3, public FlagsComponent<Algorithm> {};

template <typename Algorithm, typename T1, typename T2, typename T3, typename T4>
struct Result4 : public T1, public T2, public T3, public T4, public FlagsComponent<Algorithm> {};

/** @} */ // end of measBaseResultTemplates group

/** @} */ // end of measBaseResults group

}}} // lsst::meas::base

#endif // !LSST_MEAS_BASE_Results_h_INCLUDED
