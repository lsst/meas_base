#
# LSST Data Management System
# Copyright 2008-2017 Aura/LSST
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

# -*- lsst-python -*-

import math
import unittest

import numpy as np

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.meas.base as measBase
import lsst.utils.tests

try:
    display
    import lsst.afw.display.ds9 as ds9
except NameError:
    display = False
    displayCoeffs = False


def plantSources(bbox, kwid, sky, coordList, addPoissonNoise=True):
    """Make an exposure with stars (modelled as Gaussians)."""
    """
    @param bbox: parent bbox of exposure
    @param kwid: kernel width (and height; kernel is square)
    @param sky: amount of sky background (counts)
    @param coordList: a list of [x, y, counts, sigma], where:
        * x,y are relative to exposure origin
        * counts is the integrated counts for the star
        * sigma is the Gaussian sigma in pixels
    @param addPoissonNoise: add Poisson noise to the exposure?
    """
    # make an image with sources
    img = afwImage.ImageD(bbox)
    meanSigma = 0.0
    for coord in coordList:
        x, y, counts, sigma = coord
        meanSigma += sigma

        # make a single gaussian psf
        psf = afwDetection.GaussianPsf(kwid, kwid, sigma)

        # make an image of it and scale to the desired number of counts
        thisPsfImg = psf.computeImage(lsst.geom.PointD(int(x), int(y)))
        thisPsfImg *= counts

        # bbox a window in our image and add the fake star image
        imgSeg = img.Factory(img, thisPsfImg.getBBox())
        imgSeg += thisPsfImg
    meanSigma /= len(coordList)

    img += sky

    # add Poisson noise
    if (addPoissonNoise):
        np.random.seed(seed=1)  # make results reproducible
        imgArr = img.getArray()
        imgArr[:] = np.random.poisson(imgArr)

    # bundle into a maskedimage and an exposure
    mask = afwImage.Mask(bbox)
    var = img.convertFloat()
    img -= sky
    mimg = afwImage.MaskedImageF(img.convertFloat(), mask, var)
    exposure = afwImage.makeExposure(mimg)

    # insert an approximate psf
    psf = afwDetection.GaussianPsf(kwid, kwid, meanSigma)
    exposure.setPsf(psf)

    return exposure


class SincPhotSums(lsst.utils.tests.TestCase):

    def setUp(self):
        self.nx = 64
        self.ny = 64
        self.kwid = 15
        self.sky = 100.0
        self.val = 10000.0
        self.sigma = 4.0
        coordList = [[self.nx/2, self.ny/2, self.val, self.sigma]]

        # exposure with gaussian
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(self.nx, self.ny))
        self.expGaussPsf = plantSources(bbox, self.kwid, self.sky, coordList, addPoissonNoise=False)

        # just plain sky (ie. a constant)
        self.mimg = afwImage.MaskedImageF(lsst.geom.ExtentI(self.nx, self.ny))
        self.mimg.set(self.sky, 0x0, self.sky)
        self.expSky = afwImage.makeExposure(self.mimg)

        if display > 1:
            ds9.mtv(self.expGaussPsf)

    def tearDown(self):
        del self.mimg
        del self.expGaussPsf
        del self.expSky

    def testEllipticalGaussian(self):
        """Test measuring elliptical aperture mags for an elliptical Gaussian"""

        width, height = 200, 200
        xcen, ycen = 0.5*width, 0.5*height
        #
        # Make the object
        #
        gal = afwImage.ImageF(lsst.geom.ExtentI(width, height))
        a, b, theta = float(10), float(5), 20
        instFlux = 1e4
        I0 = instFlux/(2*math.pi*a*b)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        for y in range(height):
            for x in range(width):
                dx, dy = x - xcen, y - ycen
                u = c*dx + s*dy
                v = -s*dx + c*dy
                val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                if val < 0:
                    val = 0
                gal[x, y, afwImage.LOCAL] = val

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        del gal

        if display:
            frame = 0
            ds9.mtv(objImg, frame=frame, title="Elliptical")

        self.assertAlmostEqual(1.0, afwMath.makeStatistics(objImg.getMaskedImage().getImage(),
                                                           afwMath.SUM).getValue()/instFlux)
        #
        # Now measure some annuli
        #
        for r1, r2 in [(0., 0.45*a),
                       (0.45*a, 1.0*a),
                       (1.0*a, 2.0*a),
                       (2.0*a, 3.0*a),
                       (3.0*a, 5.0*a),
                       (3.0*a, 10.0*a),
                       ]:
            if display:                 # draw the inner and outer boundaries of the aperture
                Mxx = 1
                Myy = (b/a)**2

                mxx, mxy, myy = c**2*Mxx + s**2*Myy, c*s*(Mxx - Myy), s**2*Mxx + c**2*Myy
                for r in (r1, r2):
                    ds9.dot("@:%g,%g,%g" % (r**2*mxx, r**2*mxy, r**2*myy), xcen, ycen, frame=frame)

            center = lsst.geom.Point2D(xcen, ycen)

            # this tests tests a sync algorithm with an inner and outer radius
            # since that is no longer available from the ApertureFluxAlgorithm,
            # we will calculate the two and subtract.

            axes = afwGeom.ellipses.Axes(r2, r2*(1-b/a), math.radians(theta))
            ellipse = afwGeom.Ellipse(axes, center)
            result2 = measBase.ApertureFluxAlgorithm.computeSincFlux(objImg.getMaskedImage(), ellipse)

            axes = afwGeom.ellipses.Axes(r1, r1*(1-b/a), math.radians(theta))
            ellipse = afwGeom.Ellipse(axes, center)
            result1 = measBase.ApertureFluxAlgorithm.computeSincFlux(objImg.getMaskedImage(), ellipse)

            self.assertAlmostEqual(math.exp(-0.5*(r1/a)**2) - math.exp(-0.5*(r2/a)**2),
                                   (result2.instFlux-result1.instFlux)/instFlux, 4)


class SincCoeffTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.ellipse = afwEll.Axes(10.0, 5.0, 0.12345)
        self.radius1 = 0.1234
        self.radius2 = 4.3210
        self.inner = self.radius1/self.radius2

    def tearDown(self):
        del self.ellipse

    def assertCached(self, coeff1, coeff2):
        np.testing.assert_array_equal(coeff1.getArray(), coeff2.getArray())
        self.assertEqual(coeff1.getId(), coeff2.getId())

    def assertNotCached(self, coeff1, coeff2):
        np.testing.assert_array_equal(coeff1.getArray(), coeff2.getArray())
        self.assertNotEqual(coeff1.getId(), coeff2.getId())

    def getCoeffCircle(self, radius2):
        circle = afwEll.Axes(radius2, radius2, 0.0)
        inner = self.radius1/radius2
        coeff1 = measBase.SincCoeffsF.get(circle, inner)
        coeff2 = measBase.SincCoeffsF.get(circle, inner)
        return coeff1, coeff2

    def testNoCachingElliptical(self):
        coeff1 = measBase.SincCoeffsF.get(self.ellipse, self.inner)
        coeff2 = measBase.SincCoeffsF.get(self.ellipse, self.inner)
        self.assertNotCached(coeff1, coeff2)

    def testNoCachingCircular(self):
        coeff1, coeff2 = self.getCoeffCircle(2*self.radius2)  # not self.radius2 because that may be cached
        self.assertNotCached(coeff1, coeff2)

    def testWithCaching(self):
        measBase.SincCoeffsF.cache(self.radius1, self.radius2)
        coeff1, coeff2 = self.getCoeffCircle(self.radius2)
        self.assertCached(coeff1, coeff2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
