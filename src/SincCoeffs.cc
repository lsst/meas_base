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

#include <complex>

#include "boost/math/special_functions/bessel.hpp"
#include "boost/shared_array.hpp"

#ifndef LSST_DISABLE_SINC_PHOTOMETRY
#include "fftw3.h"
#endif

#include "lsst/meas/base/SincCoeffs.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Extent.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/math/Integrate.h"

namespace lsst {
namespace meas {
namespace base {

#ifndef LSST_DISABLE_SINC_PHOTOMETRY

namespace {

// Convenient wrapper for a Bessel function
inline double J1(double const x) { return boost::math::cyl_bessel_j(1, x); }

// sinc function
template <typename T>
inline T sinc(T const x) {
    return (x != 0.0) ? (std::sin(x) / x) : 1.0;
}

/*
 * Define a circular aperture function object g_i, cos-tapered?
 */
template <typename CoordT>
class CircularAperture {
public:
    CircularAperture(
            CoordT const radius1,    ///< inner radius of the aperture
            CoordT const radius2,    ///< outer radius of the aperture
            CoordT const taperwidth  ///< width to cosine taper from 1.0 to 0.0 (ie. 0.5*cosine period)
            )
            : _radius1(radius1),
              _radius2(radius2),
              _taperwidth1(taperwidth),
              _taperwidth2(taperwidth),
              _k1(1.0 / (2.0 * taperwidth)),
              _k2(1.0 / (2.0 * taperwidth)),
              _taperLo1(radius1 - 0.5 * taperwidth),
              _taperHi1(radius1 + 0.5 * taperwidth),
              _taperLo2(radius2 - 0.5 * taperwidth),
              _taperHi2(radius2 + 0.5 * taperwidth) {
        // if we're asked for a radius smaller than our taperwidth,
        // adjust the taper width smaller so it fits exactly
        // with smooth derivative=0 at r=0

        if (_radius1 > _radius2) {
            throw LSST_EXCEPT(
                    pex::exceptions::InvalidParameterError,
                    (boost::format("rad2 less than rad1: (rad1=%.2f, rad2=%.2f) ") % _radius1 % _radius2)
                            .str());
        }
        if (_radius1 < 0.0 || _radius2 < 0.0) {
            throw LSST_EXCEPT(
                    pex::exceptions::InvalidParameterError,
                    (boost::format("radii must be >= 0 (rad1=%.2f, rad2=%.2f) ") % _radius1 % _radius2)
                            .str());
        }

        if (_radius1 == 0) {
            _taperwidth1 = 0.0;
            _k1 = 0.0;
        }

        // if we don't have room to taper at r=0
        if (_radius1 < 0.5 * _taperwidth1) {
            _taperwidth1 = 2.0 * _radius1;
            _k1 = 1.0 / (2.0 * _taperwidth1);
        }

        // if we don't have room to taper between r1 and r2
        if ((_radius2 - _radius1) < 0.5 * (_taperwidth1 + _taperwidth2)) {
            // if we *really* don't have room ... taper1 by itself is too big
            // - set taper1,2 to be equal and split the r2-r1 range
            if ((_radius2 - _radius2) < 0.5 * _taperwidth1) {
                _taperwidth1 = _taperwidth2 = 0.5 * (_radius2 - _radius1);
                _k1 = _k2 = 1.0 / (2.0 * _taperwidth1);

                // if there's room for taper1, but not taper1 and 2
            } else {
                _taperwidth2 = _radius2 - _radius1 - _taperwidth1;
                _k2 = 1.0 / (2.0 * _taperwidth2);
            }

            _taperLo1 = _radius1 - 0.5 * _taperwidth1;
            _taperHi1 = _radius1 + 0.5 * _taperwidth1;
            _taperLo2 = _radius2 - 0.5 * _taperwidth2;
            _taperHi2 = _radius2 + 0.5 * _taperwidth2;
        }
    }

    // When called, return the throughput at the requested x,y
    // todo: replace the sinusoid taper with a band-limited
    CoordT operator()(CoordT const x, CoordT const y) const {
        CoordT const xyrad = std::sqrt(x * x + y * y);
        if (xyrad < _taperLo1) {
            return 0.0;
        } else if (xyrad >= _taperLo1 && xyrad <= _taperHi1) {
            return 0.5 * (1.0 + std::cos((geom::TWOPI * _k1) * (xyrad - _taperHi1)));
        } else if (xyrad > _taperHi1 && xyrad <= _taperLo2) {
            return 1.0;
        } else if (xyrad > _taperLo2 && xyrad <= _taperHi2) {
            return 0.5 * (1.0 + std::cos((geom::TWOPI * _k2) * (xyrad - _taperLo2)));
        } else {
            return 0.0;
        }
    }

    CoordT getRadius1() { return _radius1; }
    CoordT getRadius2() { return _radius2; }

private:
    CoordT _radius1, _radius2;
    CoordT _taperwidth1, _taperwidth2;
    CoordT _k1, _k2;  // the angular wavenumber corresponding to a cosine with wavelength 2*taperwidth
    CoordT _taperLo1, _taperHi1;
    CoordT _taperLo2, _taperHi2;
};

template <typename CoordT>
class CircApPolar : public std::unary_function<CoordT, CoordT> {
public:
    CircApPolar(double radius, double taperwidth) : _ap(CircularAperture<CoordT>(0.0, radius, taperwidth)) {}
    CoordT operator()(double r) const { return r * _ap(r, 0.0); }

private:
    CircularAperture<CoordT> _ap;
};

/*
 * Define a Sinc functor to be integrated over for Sinc interpolation
 */
template <typename IntegrandT>
class SincAperture : public std::binary_function<IntegrandT, IntegrandT, IntegrandT> {
public:
    SincAperture(CircularAperture<IntegrandT> const& ap,
                 int const ix,  // sinc center x
                 int const iy   // sinc center y
                 )
            : _ap(ap), _ix(ix), _iy(iy) {}

    IntegrandT operator()(IntegrandT const x, IntegrandT const y) const {
        double const fourierConvention = geom::PI;
        double const dx = fourierConvention * (x - _ix);
        double const dy = fourierConvention * (y - _iy);
        double const fx = sinc<double>(dx);
        double const fy = sinc<double>(dy);
        return (1.0 + _ap(x, y) * fx * fy);
    }

private:
    CircularAperture<IntegrandT> const& _ap;
    double _ix, _iy;
};

class GaussPowerFunctor : public std::binary_function<double, double, double> {
public:
    GaussPowerFunctor(double sigma) : _sigma(sigma) {}

    double operator()(double kx, double ky) const {
        double const k = ::sqrt(kx * kx + ky * ky);
        double const gauss = ::exp(-0.5 * k * k * _sigma * _sigma);
        return gauss * gauss;
    }

private:
    double _sigma;
};

std::pair<double, double> computeGaussLeakage(double const sigma) {
    GaussPowerFunctor gaussPower(sigma);

    double lim = geom::PI;

    // total power: integrate GaussPowerFunctor -inf<x<inf, -inf<y<inf (can be done analytically)
    double powerInf = geom::PI / (sigma * sigma);

    // true power: integrate GaussPowerFunctor -lim<x<lim, -lim<y<lim (must be done numerically)
    double truePower = afw::math::integrate2d(gaussPower, -lim, lim, -lim, lim, 1.0e-8);
    double trueLeak = (powerInf - truePower) / powerInf;

    // estimated power: function is circular, but coords are cartesian
    // - true power does the actual integral numerically, but we can estimate it by integrating
    //   in polar coords over lim <= radius < infinity.  The integral is analytic.
    double estLeak = ::exp(-sigma * sigma * geom::PI * geom::PI) / powerInf;

    return std::pair<double, double>(trueLeak, estLeak);
}

template <typename PixelT>
std::shared_ptr<afw::image::Image<PixelT>> calcImageRealSpace(double const rad1, double const rad2,
                                                              double const taperwidth) {
    PixelT initweight = 0.0;  // initialize the coeff values

    int log2 = static_cast<int>(::ceil(::log10(2.0 * rad2) / log10(2.0)));
    if (log2 < 3) {
        log2 = 3;
    }
    int hwid = pow(2, log2);
    int width = 2 * hwid - 1;

    int const xwidth = width;
    int const ywidth = width;

    int const x0 = -xwidth / 2;
    int const y0 = -ywidth / 2;

    // create an image to hold the coefficient image
    auto coeffImage = std::make_shared<afw::image::Image<PixelT>>(geom::ExtentI(xwidth, ywidth), initweight);
    coeffImage->setXY0(x0, y0);

    // create the aperture function object
    // determine the radius to use that makes 'radius' the effective radius of the aperture
    double tolerance = 1.0e-12;
    double dr = 1.0e-6;
    double err = 2.0 * tolerance;
    double apEff = geom::PI * rad2 * rad2;
    double radIn = rad2;
    int maxIt = 20;
    int i = 0;
    while (err > tolerance && i < maxIt) {
        CircApPolar<double> apPolar1(radIn, taperwidth);
        CircApPolar<double> apPolar2(radIn + dr, taperwidth);
        double a1 = geom::TWOPI * afw::math::integrate(apPolar1, 0.0, radIn + taperwidth, tolerance);
        double a2 = geom::TWOPI * afw::math::integrate(apPolar2, 0.0, radIn + dr + taperwidth, tolerance);
        double dadr = (a2 - a1) / dr;
        double radNew = radIn - (a1 - apEff) / dadr;
        err = (a1 - apEff) / apEff;
        radIn = radNew;
        i++;
    }
    CircularAperture<double> ap(rad1, rad2, taperwidth);

    /* ******************************* */
    // integrate over the aperture

    // the limits of the integration over the sinc aperture
    double const limit = rad2 + taperwidth;
    double const x1 = -limit;
    double const x2 = limit;
    double const y1 = -limit;
    double const y2 = limit;

    for (int iY = y0; iY != y0 + coeffImage->getHeight(); ++iY) {
        int iX = x0;
        typename afw::image::Image<PixelT>::x_iterator end = coeffImage->row_end(iY - y0);
        for (typename afw::image::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY - y0); ptr != end;
             ++ptr) {
            // create a sinc function in the CircularAperture at our location
            SincAperture<double> sincAp(ap, iX, iY);

            // integrate the sinc
            PixelT integral = afw::math::integrate2d(sincAp, x1, x2, y1, y2, 1.0e-8);

            // we actually integrated function+1.0 and now must subtract the excess volume
            // - just force it to zero in the corners
            double const dx = iX;
            double const dy = iY;
            *ptr = (std::sqrt(dx * dx + dy * dy) < xwidth / 2) ? integral - (x2 - x1) * (y2 - y1) : 0.0;
            ++iX;
        }
    }

    double sum = 0.0;
    for (int iY = y0; iY != y0 + coeffImage->getHeight(); ++iY) {
        typename afw::image::Image<PixelT>::x_iterator end = coeffImage->row_end(iY - y0);
        for (typename afw::image::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY - y0); ptr != end;
             ++ptr) {
            sum += *ptr;
        }
    }

#if 0  // debugging
    coeffImage->writeFits("cimage.fits");
#endif

    return coeffImage;
}

class FftShifter {
public:
    FftShifter(int xwid) : _xwid(xwid) {}
    int shift(int x) {
        if (x >= _xwid / 2) {
            return x - _xwid / 2;
        } else {
            return x + _xwid / 2 + 1;
        }
    }

private:
    int _xwid;
};

std::pair<double, double> rotate(double x, double y, double angle) {
    double c = ::cos(angle);
    double s = ::sin(angle);
    return std::pair<double, double>(x * c + y * s, -x * s + y * c);
}

/*  todo
 * - try sub pixel shift if it doesn't break even symmetry
 * - put values directly in an Image
 * - precompute the plan
 */

template <typename PixelT>
std::shared_ptr<afw::image::Image<PixelT>> calcImageKSpaceCplx(double const rad1, double const rad2,
                                                               double const posAng,
                                                               double const ellipticity) {
    // we only need a half-width due to symmetry
    // make the hwid 2*rad2 so we have some buffer space and round up to the next power of 2
    int log2 = static_cast<int>(::ceil(::log10(2.0 * rad2) / log10(2.0)));
    if (log2 < 3) {
        log2 = 3;
    }
    int hwid = pow(2, log2);
    int wid = 2 * hwid - 1;
    int xcen = wid / 2, ycen = wid / 2;
    FftShifter fftshift(wid);

    boost::shared_array<std::complex<double>> cimg(new std::complex<double>[wid * wid]);
    std::complex<double>* c = cimg.get();
    // fftplan args: nx, ny, *in, *out, direction, flags
    // - done in-situ if *in == *out
    fftw_plan plan = fftw_plan_dft_2d(wid, wid, reinterpret_cast<fftw_complex*>(c),
                                      reinterpret_cast<fftw_complex*>(c), FFTW_BACKWARD, FFTW_ESTIMATE);

    // compute the k-space values and put them in the cimg array
    double const twoPiRad1 = geom::TWOPI * rad1;
    double const twoPiRad2 = geom::TWOPI * rad2;
    double const scale = (1.0 - ellipticity);
    for (int iY = 0; iY < wid; ++iY) {
        int const fY = fftshift.shift(iY);
        double const ky = (static_cast<double>(iY) - ycen) / wid;

        for (int iX = 0; iX < wid; ++iX) {
            int const fX = fftshift.shift(iX);
            double const kx = static_cast<double>(iX - xcen) / wid;

            // rotate
            std::pair<double, double> coo = rotate(kx, ky, posAng);
            double kxr = coo.first;
            double kyr = coo.second;
            // rescale
            double const k = ::sqrt(kxr * kxr + scale * scale * kyr * kyr);

            double const airy1 = (rad1 > 0 ? rad1 * J1(twoPiRad1 * k) : 0.0) / k;
            double const airy2 = rad2 * J1(twoPiRad2 * k) / k;
            double const airy = airy2 - airy1;

            c[fY * wid + fX] = std::complex<double>(scale * airy, 0.0);
        }
    }
    c[0] = scale * geom::PI * (rad2 * rad2 - rad1 * rad1);

    // perform the fft and clean up after ourselves
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // put the coefficients into an image
    auto coeffImage = std::make_shared<afw::image::Image<PixelT>>(geom::ExtentI(wid, wid), 0.0);

    for (int iY = 0; iY != coeffImage->getHeight(); ++iY) {
        int iX = 0;
        typename afw::image::Image<PixelT>::x_iterator end = coeffImage->row_end(iY);
        for (typename afw::image::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY); ptr != end;
             ++ptr) {
            int fX = fftshift.shift(iX);
            int fY = fftshift.shift(iY);
            *ptr = static_cast<PixelT>(c[fY * wid + fX].real() / (wid * wid));
            iX++;
        }
    }

    // reset the origin to be the middle of the image
    coeffImage->setXY0(-wid / 2, -wid / 2);
    return coeffImage;
}

// I'm not sure why this doesn't work with DCT-I (REDFT00), the DCT should take advantage of symmetry
// and be much faster than the DFT.  It runs but the numbers are slightly off ...
// but I have to do it as real-to-halfcplx (R2HC) to get the correct numbers.

template <typename PixelT>
std::shared_ptr<afw::image::Image<PixelT>> calcImageKSpaceReal(double const rad1, double const rad2) {
    // we only need a half-width due to symmertry
    // make the hwid 2*rad2 so we have some buffer space and round up to the next power of 2
    int log2 = static_cast<int>(::ceil(::log10(2.0 * rad2) / log10(2.0)));
    if (log2 < 3) {
        log2 = 3;
    }
    int hwid = pow(2, log2);
    int wid = 2 * hwid - 1;
    int xcen = wid / 2, ycen = wid / 2;
    FftShifter fftshift(wid);

    boost::shared_array<double> cimg(new double[wid * wid]);
    double* c = cimg.get();
    // fftplan args: nx, ny, *in, *out, kindx, kindy, flags
    // - done in-situ if *in == *out
    fftw_plan plan = fftw_plan_r2r_2d(wid, wid, c, c, FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);

    // compute the k-space values and put them in the cimg array
    double const twoPiRad1 = geom::TWOPI * rad1;
    double const twoPiRad2 = geom::TWOPI * rad2;
    for (int iY = 0; iY < wid; ++iY) {
        int const fY = fftshift.shift(iY);
        double const ky = (static_cast<double>(iY) - ycen) / wid;

        for (int iX = 0; iX < wid; ++iX) {
            int const fX = fftshift.shift(iX);

            // emacs indent breaks if this isn't separte
            double const iXcen = static_cast<double>(iX - xcen);
            double const kx = iXcen / wid;

            double const k = ::sqrt(kx * kx + ky * ky);
            double const airy1 = (rad1 > 0 ? rad1 * J1(twoPiRad1 * k) : 0.0) / k;
            double const airy2 = rad2 * J1(twoPiRad2 * k) / k;
            double const airy = airy2 - airy1;
            c[fY * wid + fX] = airy;
        }
    }
    int fxy = fftshift.shift(wid / 2);
    c[fxy * wid + fxy] = geom::PI * (rad2 * rad2 - rad1 * rad1);

    // perform the fft and clean up after ourselves
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // put the coefficients into an image
    auto coeffImage = std::make_shared<afw::image::Image<PixelT>>(geom::ExtentI(wid, wid), 0.0);

    for (int iY = 0; iY != coeffImage->getHeight(); ++iY) {
        int iX = 0;
        typename afw::image::Image<PixelT>::x_iterator end = coeffImage->row_end(iY);
        for (typename afw::image::Image<PixelT>::x_iterator ptr = coeffImage->row_begin(iY); ptr != end;
             ++ptr) {
            // now need to reflect the quadrant we solved to the other three
            int fX = iX < hwid ? hwid - iX - 1 : iX - hwid + 1;
            int fY = iY < hwid ? hwid - iY - 1 : iY - hwid + 1;
            *ptr = static_cast<PixelT>(c[fY * wid + fX] / (wid * wid));
            iX++;
        }
    }

    // reset the origin to be the middle of the image
    coeffImage->setXY0(-wid / 2, -wid / 2);
    return coeffImage;
}

}  // namespace

#endif // !LSST_DISABLE_SINC_PHOTOMETRY

template <typename PixelT>
SincCoeffs<PixelT>& SincCoeffs<PixelT>::getInstance() {
    static SincCoeffs<PixelT> instance;
    return instance;
}

template <typename PixelT>
void SincCoeffs<PixelT>::cache(float r1, float r2) {
#ifdef LSST_DISABLE_SINC_PHOTOMETRY
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Sinc photometry has been disabled at compile-time.");
#else
    if (r1 < 0.0 || r2 < r1) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          (boost::format("Invalid r1,r2 = %f,%f") % r1 % r2).str());
    }
    double const innerFactor = r1 / r2;
    afw::geom::ellipses::Axes axes(r2, r2, 0.0);
    if (!getInstance()._lookup(axes, innerFactor)) {
        PTR(typename SincCoeffs<PixelT>::CoeffT) coeff = calculate(axes, innerFactor);
        coeff->markPersistent();
        getInstance()._cache[r2][innerFactor] = coeff;
    }
#endif
}

template <typename PixelT>
CONST_PTR(typename SincCoeffs<PixelT>::CoeffT)
SincCoeffs<PixelT>::get(afw::geom::ellipses::Axes const& axes, float const innerFactor) {
#ifdef LSST_DISABLE_SINC_PHOTOMETRY
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Sinc photometry has been disabled at compile-time.");
#else
    CONST_PTR(CoeffT) coeff = getInstance()._lookup(axes, innerFactor);
    return coeff ? coeff : calculate(axes, innerFactor);
#endif
}

template <typename PixelT>
CONST_PTR(typename SincCoeffs<PixelT>::CoeffT)
SincCoeffs<PixelT>::_lookup(afw::geom::ellipses::Axes const& axes, double const innerFactor) const {
#ifdef LSST_DISABLE_SINC_PHOTOMETRY
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Sinc photometry has been disabled at compile-time.");
#else
    if (innerFactor < 0.0 || innerFactor > 1.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          (boost::format("innerFactor = %f is not between 0 and 1") % innerFactor).str());
    }

    CONST_PTR(typename SincCoeffs<PixelT>::CoeffT) const null = CONST_PTR(SincCoeffs<PixelT>::CoeffT)();

    // We only cache circular apertures
    if (!FuzzyCompare<float>().isEqual(axes.getA(), axes.getB())) {
        return null;
    }
    typename CoeffMapMap::const_iterator iter1 = _cache.find(axes.getA());
    if (iter1 == _cache.end()) {
        return null;
    }
    typename CoeffMap::const_iterator iter2 = iter1->second.find(innerFactor);
    return (iter2 == iter1->second.end()) ? null : iter2->second;
#endif
}

template <typename PixelT>
PTR(typename SincCoeffs<PixelT>::CoeffT)
SincCoeffs<PixelT>::calculate(afw::geom::ellipses::Axes const& axes, double const innerFactor) {
#ifdef LSST_DISABLE_SINC_PHOTOMETRY
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Sinc photometry has been disabled at compile-time.");
#else
    if (innerFactor < 0.0 || innerFactor > 1.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          (boost::format("innerFactor = %f is not between 0 and 1") % innerFactor).str());
    }

    // Kspace-real is fastest, but only slightly faster than kspace cplx
    // but real won't work for elliptical apertures due to symmetries assumed for real transform

    double const rad1 = axes.getA() * innerFactor;
    double const rad2 = axes.getA();
    // if there's no angle and no ellipticity
    if (FuzzyCompare<float>().isEqual(axes.getA(), axes.getB())) {
        // here we call the real transform
        return calcImageKSpaceReal<PixelT>(rad1, rad2);
    } else {
        // here we call the complex transform
        double const ellipticity = 1.0 - axes.getB() / axes.getA();
        return calcImageKSpaceCplx<PixelT>(rad1, rad2, axes.getTheta(), ellipticity);
    }
#endif
}

template class SincCoeffs<float>;
template class SincCoeffs<double>;

}  // namespace base
}  // namespace meas
}  // namespace lsst
