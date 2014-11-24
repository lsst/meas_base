#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

import unittest
import numpy

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
from lsst.meas.base import ApertureFluxAlgorithm
import lsst.meas.base.tests

numpy.random.seed(1234567)

# n.b. Some tests here depend on the noise realization in the test data
# or from the numpy random number generator.
# For the current test data and seed value, they pass, but they may not
# if the test data is regenerated or the seed value changes.  I've marked
# these with an "rng dependent" comment.  In most cases, they test that
# the measured flux lies within 2 sigma of the correct value, which we
# should expect to fail sometimes.

class ApertureFluxTestCase(lsst.utils.tests.TestCase):
    """Test case for the ApertureFlux algorithm base class
    """

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(20, -100), lsst.afw.geom.Point2I(100, -20))
        self.exposure = lsst.afw.image.ExposureF(self.bbox)
        self.exposure.getMaskedImage().getImage().set(1.0)
        self.exposure.getMaskedImage().getVariance().set(0.25)
        self.ctrl = ApertureFluxAlgorithm.Control()

    def tearDown(self):
        del self.bbox
        del self.exposure

    def computeNaiveArea(self, position, radius):
        """Return the area of a circular aperture with the given position and radius, according to
        the 'naive' definition of the aperture - just test whether the center of each pixel is within
        the circle.
        """
        x, y = numpy.meshgrid(numpy.arange(self.bbox.getBeginX(), self.bbox.getEndX()),
                              numpy.arange(self.bbox.getBeginY(), self.bbox.getEndY()))
        return ((x - position.getX())**2 + (y - position.getY())**2 <= radius**2).sum()

    def testNaive(self):
        positions = [lsst.afw.geom.Point2D(60.0, -60.0),
                     lsst.afw.geom.Point2D(60.5, -60.0),
                     lsst.afw.geom.Point2D(60.0, -60.5),
                     lsst.afw.geom.Point2D(60.5, -60.5)]
        radii = [12.0, 17.0]
        for position in positions:
            for radius in radii:
                ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(radius, radius, 0.0),
                                                         position)
                area = self.computeNaiveArea(position, radius)
                # test that this isn't the same as the sinc flux
                self.assertNotClose(
                    ApertureFluxAlgorithm.computeSincFlux(self.exposure.getMaskedImage().getImage(),
                                                          ellipse, self.ctrl).flux,
                    area
                )
                # test that all the ways we could invoke naive flux measurement produce the expected result
                def check(method, image):
                    result = method(image, ellipse, self.ctrl)
                    self.assertClose(result.flux, area)
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED))
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED))
                    if hasattr(image, "getVariance"):
                        self.assertClose(result.fluxSigma, (area*0.25)**0.5)
                    else:
                        self.assertTrue(numpy.isnan(result.fluxSigma))
                check(ApertureFluxAlgorithm.computeNaiveFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeNaiveFlux, self.exposure.getMaskedImage().getImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage().getImage())
        # test failure conditions when the aperture itself is truncated
        invalid = ApertureFluxAlgorithm.computeNaiveFlux(
            self.exposure.getMaskedImage().getImage(),
            lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(12.0, 12.0),
                                           lsst.afw.geom.Point2D(25.0, -60.0)),
            self.ctrl
            )
        self.assertTrue(invalid.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED))
        self.assertFalse(invalid.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED))
        self.assertTrue(numpy.isnan(invalid.flux))

    def testSinc(self):
        positions = [lsst.afw.geom.Point2D(60.0, -60.0),
                     lsst.afw.geom.Point2D(60.5, -60.0),
                     lsst.afw.geom.Point2D(60.0, -60.5),
                     lsst.afw.geom.Point2D(60.5, -60.5)]
        radii = [7.0, 9.0]
        for position in positions:
            for radius in radii:
                ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(radius, radius, 0.0),
                                                         position)
                area = ellipse.getCore().getArea()
                # test that this isn't the same as the naive flux
                self.assertNotClose(
                    ApertureFluxAlgorithm.computeNaiveFlux(self.exposure.getMaskedImage().getImage(),
                                                           ellipse, self.ctrl).flux,
                    area
                )
                # test that all the ways we could invoke sinc flux measurement produce the expected result
                def check(method, image):
                    result = method(image, ellipse, self.ctrl)
                    self.assertClose(result.flux, area, rtol=1E-3)
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED))
                    self.assertFalse(result.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED))
                    if hasattr(image, "getVariance"):
                        self.assertFalse(numpy.isnan(result.fluxSigma))
                    else:
                        self.assertTrue(numpy.isnan(result.fluxSigma))
                check(ApertureFluxAlgorithm.computeSincFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeSincFlux, self.exposure.getMaskedImage().getImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage())
                check(ApertureFluxAlgorithm.computeFlux, self.exposure.getMaskedImage().getImage())
        # test failure conditions when the aperture itself is truncated
        invalid1 = ApertureFluxAlgorithm.computeSincFlux(
            self.exposure.getMaskedImage().getImage(),
            lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(9.0, 9.0),
                                           lsst.afw.geom.Point2D(25.0, -60.0)),
            self.ctrl
            )
        self.assertTrue(invalid1.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED))
        self.assertTrue(invalid1.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED))
        self.assertTrue(numpy.isnan(invalid1.flux))
        # test failure conditions when the aperture is not truncated, but the sinc coeffs are
        invalid2 = ApertureFluxAlgorithm.computeSincFlux(
            self.exposure.getMaskedImage().getImage(),
            lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(9.0, 9.0),
                                           lsst.afw.geom.Point2D(30.0, -60.0)),
            self.ctrl
            )
        self.assertFalse(invalid2.getFlag(ApertureFluxAlgorithm.APERTURE_TRUNCATED))
        self.assertTrue(invalid2.getFlag(ApertureFluxAlgorithm.SINC_COEFFS_TRUNCATED))
        self.assertFalse(numpy.isnan(invalid2.flux))

class CircularApertureFluxTestCase(lsst.meas.base.tests.AlgorithmTestCase):
    """Test case for the CircularApertureFlux algorithm/plugin
    """

    def setUp(self):
        lsst.meas.base.tests.AlgorithmTestCase.setUp(self)

    def tearDown(self):
        lsst.meas.base.tests.AlgorithmTestCase.tearDown(self)

    def testSingleFramePlugin(self):
        config = lsst.meas.base.SingleFrameMeasurementConfig()
        config.plugins["base_CircularApertureFlux"].maxSincRadius = 20
        ctrl = config.plugins["base_CircularApertureFlux"].makeControl()
        measCat = self.runSingleFrameMeasurementTask("base_CircularApertureFlux", config=config)
        metadata = measCat.getMetadata()
        radii = metadata.get("base_CircularApertureFlux_radii")
        #self.assertEqual(list(radii), list(ctrl.radii))
        for record in measCat:
            lastFlux = 0.0
            lastFluxSigma = 0.0
            for n, radius in enumerate(radii):
                # Test that the flags are what we expect
                if radius <= ctrl.maxSincRadius:
                    self.assertFalse(record.get("base_CircularApertureFlux_flag_%d" % n))
                    self.assertFalse(record.get("base_CircularApertureFlux_flag_apertureTruncated_%d" % n))
                    self.assertEqual(
                        record.get("base_CircularApertureFlux_flag_sincCoeffsTruncated_%d" % n),
                        radius > 12
                    )
                else:
                    self.assertFalse(record.get("base_CircularApertureFlux_flag_sincCoeffsTruncated_%d" % n))
                    self.assertEqual(record.get("base_CircularApertureFlux_flag_%d" % n), radius > 50)
                    self.assertEqual(record.get("base_CircularApertureFlux_flag_apertureTruncated_%d" % n),
                                     radius > 50)
                # Test that the fluxes and uncertainties increase as we increase the apertures, or that
                # they match the true flux within 2 sigma.  This is just a test as to whether the values
                # are reasonable.  As to whether the values are exactly correct, we rely on the tests on
                # ApertureFluxAlgorithm's static methods, as the way the plugins code calls that is
                # extremely simple, so if the results we get are reasonable, it's hard to imagine
                # how they could be incorrect if ApertureFluxAlgorithm's tests are valid.
                currentFlux = record.get("base_CircularApertureFlux_flux_%d" % n)
                currentFluxSigma = record.get("base_CircularApertureFlux_fluxSigma_%d" % n)
                if not record.get("base_CircularApertureFlux_flag_%d" % n):
                    self.assertTrue(currentFlux > lastFlux
                                    or (record.get("truth_flux") - currentFlux) < 2*currentFluxSigma)
                    self.assertGreater(currentFluxSigma, lastFluxSigma)
                    lastFlux = currentFlux
                    lastFluxSigma = currentFluxSigma
                else:
                    self.assertTrue(numpy.isnan(currentFlux))
                    self.assertTrue(numpy.isnan(currentFluxSigma))

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ApertureFluxTestCase)
    suites += unittest.makeSuite(CircularApertureFluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
