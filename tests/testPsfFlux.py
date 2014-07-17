#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import os
import unittest
import numpy

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.base.tests

numpy.random.seed(1234567)

# n.b. Some tests here depend on the noise realization in the test data
# or from the numpy random number generator.
# For the current test data and seed value, they pass, but they may not
# if the test data is regenerated or the seed value changes.  I've marked
# these with an "rng dependent" comment.  In most cases, they test that
# the measured flux lies within 2 sigma of the correct value, which we
# should expect to fail sometimes.

class PsfFluxTestCase(lsst.meas.base.tests.AlgorithmTestCase):
    """Test case for the PsfFlux algorithm/plugins
    """

    def setUp(self):
        lsst.meas.base.tests.AlgorithmTestCase.setUp(self)
        self.record = self.truth[0]

    def tearDown(self):
        lsst.meas.base.tests.AlgorithmTestCase.tearDown(self)
        del self.record

    def testOverloads(self):
        """Verify that we get reasonable answers when we call apply() directly on an isolated source,
        and that these results are identical regardless of which overload we use.
        """
        result1 = lsst.meas.base.PsfFluxAlgorithm.apply(self.calexp, self.record.get(self.centroidKey))
        result2 = lsst.meas.base.PsfFluxAlgorithm.apply(
            self.calexp,
            lsst.meas.base.PsfFluxAlgorithm.Input(self.record.getFootprint(),
                                                  self.record.get(self.centroidKey))
            )
        self.assertEqual(result1.flux, result2.flux)
        self.assertEqual(result1.fluxSigma, result2.fluxSigma)
        self.assertFalse(result1.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_PSF))
        self.assertFalse(result1.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_GOOD_PIXELS))
        self.assertFalse(result1.getFlag(lsst.meas.base.PsfFluxAlgorithm.EDGE))
        self.assertFalse(result2.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_PSF))
        self.assertFalse(result2.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_GOOD_PIXELS))
        self.assertFalse(result2.getFlag(lsst.meas.base.PsfFluxAlgorithm.EDGE))
        # rng dependent
        self.assertClose(result1.flux, self.record.get(self.fluxKey), atol=2*result1.fluxSigma)

    def testMasking(self):
        badPoint = lsst.afw.geom.Point2I(self.record.get(self.centroidKey)) + lsst.afw.geom.Extent2I(3,4)
        imageArray = self.calexp.getMaskedImage().getImage().getArray()
        maskArray = self.calexp.getMaskedImage().getMask().getArray()
        badMask = self.calexp.getMaskedImage().getMask().getPlaneBitMask("BAD")
        imageArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] = numpy.inf
        maskArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] |= badMask
        # Should get an infinite value exception, because we didn't mask that one pixel
        self.assertRaises(
            lsst.meas.base.PixelValueError,
            lsst.meas.base.PsfFluxAlgorithm.apply,
            self.calexp,
            self.record.get(self.centroidKey)
            )
        # If we do mask it, we should get a reasonable result
        ctrl = lsst.meas.base.PsfFluxAlgorithm.Control()
        ctrl.badMaskPlanes = ["BAD"]
        result = lsst.meas.base.PsfFluxAlgorithm.apply(
            self.calexp,
            self.record.get(self.centroidKey),
            ctrl
            )
        # rng dependent
        self.assertClose(result.flux, self.record.get(self.fluxKey), atol=2*result.fluxSigma)
        # If we mask the whole image, we should get a MeasurementError
        maskArray[:,:] |= badMask
        self.assertRaises(
            lsst.meas.base.MeasurementError,
            lsst.meas.base.PsfFluxAlgorithm.apply,
            self.calexp,
            self.record.get(self.centroidKey),
            ctrl
            )

    def testSubImage(self):
        """Test that we don't get confused by images with nonzero xy0, and that the EDGE flag is set
        when it should be.
        """
        psfImage = self.calexp.getPsf().computeImage(self.record.get(self.centroidKey))
        bbox = psfImage.getBBox(lsst.afw.image.PARENT)
        bbox.grow(-1)
        subExposure = self.calexp.Factory(self.calexp, bbox)
        result = lsst.meas.base.PsfFluxAlgorithm.apply(
            subExposure,
            self.record.get(self.centroidKey),
            )
        self.assertClose(result.flux, self.record.get(self.fluxKey), atol=result.fluxSigma)
        self.assertTrue(result.getFlag(lsst.meas.base.PsfFluxAlgorithm.EDGE))

    def testNoPsf(self):
        """Test that we raise MeasurementError when there's no PSF.
        """
        self.calexp.setPsf(None)
        self.assertRaises(
            lsst.meas.base.MeasurementError,
            lsst.meas.base.PsfFluxAlgorithm.apply,
            self.calexp,
            self.record.get(self.centroidKey),
            )

    def testMonteCarlo(self):
        """Test that we get exactly the right answer on an ideal sim with no noise, and that
        the reported uncertainty agrees with a Monte Carlo test of the noise.

        We'll do this test in double precision, just to make sure that template also works.
        """
        psf = self.calexp.getPsf()
        position = self.record.get(self.centroidKey)
        image = psf.computeImage(position)
        exposure1 = lsst.afw.image.ExposureD(lsst.afw.image.MaskedImageD(image))
        exposure1.setPsf(psf)
        footprint = lsst.afw.detection.Footprint(image.getBBox(lsst.afw.image.PARENT))
        result1 = lsst.meas.base.PsfFluxAlgorithm.apply(exposure1, position)
        self.assertClose(result1.flux, 1.0)
        self.assertClose(result1.fluxSigma, 0.0)
        width = image.getWidth()
        height = image.getHeight();
        for pixelNoiseSigma in (0.001, 0.01, 0.1):
            fluxes = []
            fluxSigmas = []
            nSamples = 1000
            for repeat in xrange(nSamples):
                exposure2 = exposure1.clone()
                exposure2.getMaskedImage().getVariance().set(pixelNoiseSigma**2)
                noise = numpy.random.randn(height, width)* pixelNoiseSigma
                exposure2.getMaskedImage().getImage().getArray()[:,:] += noise
                result2 = lsst.meas.base.PsfFluxAlgorithm.apply(exposure2, position)
                fluxes.append(result2.flux)
                fluxSigmas.append(result2.fluxSigma)
            fluxMean = numpy.mean(fluxes)
            fluxSigmaMean = numpy.mean(fluxSigmas)
            fluxStandardDeviation = numpy.std(fluxes)
            self.assertClose(fluxSigmaMean, fluxStandardDeviation, rtol=0.10)   # rng dependent
            self.assertLess(fluxMean - 1.0, 2.0*fluxSigmaMean / nSamples**0.5)   # rng dependent

    def testSingleFramePlugin(self):
        measCat = self.runSingleFrameMeasurementTask("base_PsfFlux")
        for record in measCat:
            self.assertFalse(record.get("base_PsfFlux_flag"))
            self.assertFalse(record.get("base_PsfFlux_flag_noPsf"))
            self.assertFalse(record.get("base_PsfFlux_flag_noGoodPixels"))
            self.assertFalse(record.get("base_PsfFlux_flag_edge"))
            if record.get("truth_isStar"):
                self.assertClose(record.get("truth_flux"), record.get("base_PsfFlux_flux"),
                                 atol=3*record.get("base_PsfFlux_fluxSigma"))

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
