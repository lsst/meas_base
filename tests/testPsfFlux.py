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

    def prep(self, ctrl=None):
        """Construct an algorithm (finishing the schema in the process),
        create a record it can fill, and return both.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.PsfFluxControl()
        schema = self.mapper.getOutputSchema()
        algorithm = lsst.meas.base.PsfFluxAlgorithm(ctrl, "base_PsfFlux", schema)
        measCat = lsst.afw.table.SourceCatalog(schema)
        measCat.defineCentroid("truth")
        measRecord = measCat.addNew()
        measRecord.assign(self.truth[0], self.mapper)
        return algorithm, measRecord

    def setUp(self):
        lsst.meas.base.tests.AlgorithmTestCase.setUp(self)
        self.mapper = lsst.afw.table.SchemaMapper(self.truth.schema)
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.mapper.addMapping(self.truth.schema.find("truth_x").key)
        self.mapper.addMapping(self.truth.schema.find("truth_y").key)
        self.mapper.addMapping(self.truth.schema.find("truth_flux").key)

    def tearDown(self):
        lsst.meas.base.tests.AlgorithmTestCase.tearDown(self)
        del self.mapper

    def testMasking(self):
        badPoint = lsst.afw.geom.Point2I(self.truth[0].get(self.centroidKey)) + lsst.afw.geom.Extent2I(3,4)
        imageArray = self.calexp.getMaskedImage().getImage().getArray()
        maskArray = self.calexp.getMaskedImage().getMask().getArray()
        badMask = self.calexp.getMaskedImage().getMask().getPlaneBitMask("BAD")
        imageArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] = numpy.inf
        maskArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] |= badMask
        # Should get an infinite value exception, because we didn't mask that one pixel
        algorithm, measRecord = self.prep()
        self.assertRaises(lsst.meas.base.PixelValueError, algorithm.measure, measRecord, self.calexp)
        # If we do mask it, we should get a reasonable result
        ctrl = lsst.meas.base.PsfFluxControl()
        ctrl.badMaskPlanes = ["BAD"]
        algorithm, measRecord = self.prep(ctrl)
        algorithm.measure(measRecord, self.calexp)
        # rng dependent
        self.assertClose(measRecord.get("base_PsfFlux_flux"), self.truth[0].get(self.fluxKey),
                         atol=3*measRecord.get("base_PsfFlux_fluxSigma"))
        # If we mask the whole image, we should get a MeasurementError
        maskArray[:,:] |= badMask
        self.assertRaises(lsst.meas.base.MeasurementError, algorithm.measure, measRecord, self.calexp)

    def testSubImage(self):
        """Test that we don't get confused by images with nonzero xy0, and that the EDGE flag is set
        when it should be.
        """
        psfImage = self.calexp.getPsf().computeImage(self.truth[0].get(self.centroidKey))
        bbox = psfImage.getBBox()
        bbox.grow(-1)
        subExposure = self.calexp.Factory(self.calexp, bbox, lsst.afw.image.LOCAL)
        algorithm, measRecord = self.prep()
        algorithm.measure(measRecord, subExposure)
        self.assertClose(measRecord.get("base_PsfFlux_flux"), self.truth[0].get(self.fluxKey),
                         atol=3*measRecord.get("base_PsfFlux_fluxSigma"))
        self.assertTrue(measRecord.get("base_PsfFlux_flag_edge"))

    def testNoPsf(self):
        """Test that we raise MeasurementError when there's no PSF.
        """
        self.calexp.setPsf(None)
        algorithm, measRecord = self.prep()
        self.assertRaises(lsst.meas.base.FatalAlgorithmError, algorithm.measure, measRecord, self.calexp)

    def testMonteCarlo(self):
        """Test that we get exactly the right answer on an ideal sim with no noise, and that
        the reported uncertainty agrees with a Monte Carlo test of the noise.

        We'll do this test in double precision, just to make sure that template also works.
        """
        psf = self.calexp.getPsf()
        position = self.truth[0].get(self.centroidKey)
        image = psf.computeImage(position)
        exposure1 = lsst.afw.image.ExposureF(lsst.afw.image.MaskedImageF(image.convertF()))
        exposure1.setPsf(psf)
        footprint = lsst.afw.detection.Footprint(image.getBBox())
        algorithm, measRecord = self.prep()
        algorithm.measure(measRecord, exposure1)
        self.assertClose(measRecord.get("base_PsfFlux_flux"), 1.0, rtol=1E-8)
        self.assertClose(measRecord.get("base_PsfFlux_fluxSigma"), 0.0, rtol=1E-8)
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
                algorithm.measure(measRecord, exposure2)
                fluxes.append(measRecord.get("base_PsfFlux_flux"))
                fluxSigmas.append(measRecord.get("base_PsfFlux_fluxSigma"))
            fluxMean = numpy.mean(fluxes)
            fluxSigmaMean = numpy.mean(fluxSigmas)
            fluxStandardDeviation = numpy.std(fluxes)
            self.assertClose(fluxSigmaMean, fluxStandardDeviation, rtol=0.10)   # rng dependent
            self.assertLess(fluxMean - 1.0, 2.0*fluxSigmaMean / nSamples**0.5)   # rng dependent

    def testSingleFramePlugin(self):
        measCat = self.runSingleFrameMeasurementTask("base_PsfFlux")
        for record in measCat:
            self.assertFalse(record.get("base_PsfFlux_flag"))
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
