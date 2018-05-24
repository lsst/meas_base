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

#ifndef LSST_MEAS_BASE_SincCoeffs_h_INCLUDED
#define LSST_MEAS_BASE_SincCoeffs_h_INCLUDED

#include <map>

#include "lsst/afw/image/Image.h"
#include "lsst/afw/geom/ellipses/Axes.h"

namespace lsst {
namespace meas {
namespace base {

/**
 * A singleton to calculate and cache the coefficients for sinc photometry
 *
 * Caching is only performed for circular apertures (because elliptical
 * apertures are assumed to be generated dynamically, and hence not expected
 * to recur).  Caching must be explicitly requested for a particular circular
 * aperture (using the 'cache' method).
 */
template <typename PixelT>
class SincCoeffs {
public:
    typedef afw::image::Image<PixelT> CoeffT;

    /**
     * Cache the coefficients for a particular aperture
     *
     * The aperture is a circular annulus.
     */
    static void cache(float rInner, float rOuter);

    /**
     * Get the coefficients for an aperture
     *
     * Coefficients are retrieved from the cache, if available; otherwise they will be generated.
     */
    static PTR(CoeffT const)
            get(afw::geom::ellipses::Axes const& outerEllipse, float const innerRadiusFactor = 0.0);

    /// Calculate the coefficients for an aperture
    static PTR(CoeffT)
            calculate(afw::geom::ellipses::Axes const& outerEllipse, double const innerFactor = 0.0);

private:
    // A comparison function that doesn't require equality closer than machine epsilon
    template <typename T>
    struct FuzzyCompare {
        bool operator()(T x, T y) const {
            if (isEqual(x, y)) {
                return false;
            }
            return (x - y < 0) ? true : false;
        }
        bool isEqual(T x, T y) const { return ::fabs(x - y) < std::numeric_limits<T>::epsilon(); }
    };

    typedef std::map<float, PTR(CoeffT), FuzzyCompare<float> > CoeffMap;
    typedef std::map<float, CoeffMap, FuzzyCompare<float> > CoeffMapMap;
    SincCoeffs() : _cache(){};
    SincCoeffs(SincCoeffs const&);      // unimplemented: singleton
    void operator=(SincCoeffs const&);  // unimplemented: singleton

    static SincCoeffs& getInstance();

    /*
     * Search the cache for coefficients for an aperture
     *
     * If the coefficients are not cached, a null shared_ptr will be returned.
     */
    PTR(CoeffT const)
    _lookup(afw::geom::ellipses::Axes const& outerEllipse, double const innerRadiusFactor = 0.0) const;

    CoeffMapMap _cache;  //< Cache of coefficients
};

}  // namespace base
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_BASE_SincCoeffs_h_INCLUDED
