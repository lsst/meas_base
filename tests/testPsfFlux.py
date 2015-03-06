#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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

import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.utils.tests as utilsTests

from lsst.meas.base.tests import AlgorithmTestCase, TransformTestCase

class PsfFluxTestCase(AlgorithmTestCase):
    def prep(self, ctrl=None):
        """Construct an algorithm (finishing the schema in the process),
        create a record it can fill, and return both.
        """
        if ctrl is None:
            ctrl = measBase.PsfFluxControl()
        schema = self.mapper.getOutputSchema()
        schema.getAliasMap().set("slot_Centroid", "truth")
        algorithm = measBase.PsfFluxAlgorithm(ctrl, "base_PsfFlux", schema)
        measCat = afwTable.SourceCatalog(schema)
        measRecord = measCat.addNew()
        measRecord.assign(self.truth[0], self.mapper)
        return algorithm, measRecord

    def setUp(self):
        AlgorithmTestCase.setUp(self)
        self.mapper = afwTable.SchemaMapper(self.truth.schema)
        self.mapper.addMinimalSchema(afwTable.SourceTable.makeMinimalSchema())
        self.mapper.addMapping(self.truth.schema.find("truth_x").key)
        self.mapper.addMapping(self.truth.schema.find("truth_y").key)
        self.mapper.addMapping(self.truth.schema.find("truth_flux").key)

    def tearDown(self):
        AlgorithmTestCase.tearDown(self)
        del self.mapper

    def testMasking(self):
        badPoint = afwGeom.Point2I(self.truth[0].get(self.centroidKey)) + afwGeom.Extent2I(3,4)
        imageArray = self.calexp.getMaskedImage().getImage().getArray()
        maskArray = self.calexp.getMaskedImage().getMask().getArray()
        badMask = self.calexp.getMaskedImage().getMask().getPlaneBitMask("BAD")
        imageArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] = numpy.inf
        maskArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] |= badMask
        # Should get an infinite value exception, because we didn't mask that one pixel
        algorithm, measRecord = self.prep()
        self.assertRaises(measBase.PixelValueError, algorithm.measure, measRecord, self.calexp)
        # If we do mask it, we should get a reasonable result
        ctrl = measBase.PsfFluxControl()
        ctrl.badMaskPlanes = ["BAD"]
        algorithm, measRecord = self.prep(ctrl)
        algorithm.measure(measRecord, self.calexp)
        # rng dependent
        self.assertClose(measRecord.get("base_PsfFlux_flux"), self.truth[0].get(self.fluxKey),
                         atol=3*measRecord.get("base_PsfFlux_fluxSigma"))
        # If we mask the whole image, we should get a MeasurementError
        maskArray[:,:] |= badMask
        self.assertRaises(measBase.MeasurementError, algorithm.measure, measRecord, self.calexp)

    def testSubImage(self):
        """Test that we don't get confused by images with nonzero xy0, and that the EDGE flag is set
        when it should be.
        """
        psfImage = self.calexp.getPsf().computeImage(self.truth[0].get(self.centroidKey))
        bbox = psfImage.getBBox()
        bbox.grow(-1)
        subExposure = self.calexp.Factory(self.calexp, bbox, afwImage.LOCAL)
        algorithm, measRecord = self.prep()
        algorithm.measure(measRecord, subExposure)
        self.assertClose(measRecord.get("base_PsfFlux_flux"), self.truth[0].get(self.fluxKey),
                         atol=3*measRecord.get("base_PsfFlux_fluxSigma"))
        self.assertTrue(measRecord.get("base_PsfFlux_flag_edge"))

    def testNoPsf(self):
        """Test that we raise FatalAlgorithmError when there's no PSF.
        """
        self.calexp.setPsf(None)
        algorithm, measRecord = self.prep()
        self.assertRaises(measBase.FatalAlgorithmError, algorithm.measure, measRecord, self.calexp)

    def testMonteCarlo(self):
        """Test that we get exactly the right answer on an ideal sim with no noise, and that
        the reported uncertainty agrees with a Monte Carlo test of the noise.

        We'll do this test in double precision, just to make sure that template also works.
        """
        psf = self.calexp.getPsf()
        position = self.truth[0].get(self.centroidKey)
        image = psf.computeImage(position)
        exposure1 = afwImage.ExposureF(afwImage.MaskedImageF(image.convertF()))
        exposure1.setPsf(psf)
        footprint = afwDetection.Footprint(image.getBBox())
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
            self.assertLess(fluxMean - 1.0, 2.0*fluxSigmaMean / nSamples**0.5)  # rng dependent

    def testSingleFramePlugin(self):
        measCat = self.runSingleFrameMeasurementTask("base_PsfFlux")
        for record in measCat:
            self.assertFalse(record.get("base_PsfFlux_flag"))
            self.assertFalse(record.get("base_PsfFlux_flag_noGoodPixels"))
            self.assertFalse(record.get("base_PsfFlux_flag_edge"))
            if record.get("truth_isStar"):
                self.assertClose(record.get("truth_flux"), record.get("base_PsfFlux_flux"),
                                 atol=3*record.get("base_PsfFlux_fluxSigma"))


class PsfFluxTransformTestCase(TransformTestCase):
    controlClass = measBase.PsfFluxControl
    algorithmClass = measBase.PsfFluxAlgorithm
    transformClass = measBase.PsfFluxTransform

    def testTransform(self):
        flux, fluxSigma = 10, 1 # Arbitrary values for testing
        flagNames = ('_flag', '_flag_noGoodPixels', '_flag_edge') # Used by this algorithm

        for flag in (True, False):
            record = self.inputCat.addNew()
            record[self.name + '_flux'] = flux
            record[self.name + '_fluxSigma'] = fluxSigma
            for flagName in flagNames:
                record.set(self.name + flagName, flag)

        self._runTransform()

        # We are not testing the Calib itself; we assume that it is correct
        # when converting flux to magnitude, and merely check that the
        # transform has used it properly
        mag, magErr = self.calexp.getCalib().getMagnitude(flux, fluxSigma)
        for inSrc, outSrc in zip(self.inputCat, self.outputCat):
            self.assertEqual(outSrc[self.name + '_mag'], mag)
            self.assertEqual(outSrc[self.name + '_magErr'], magErr)
            for flagName in flagNames:
                self.assertEqual(outSrc.get(self.name + flagName), inSrc.get(self.name + flagName))

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PsfFluxTestCase)
    suites += unittest.makeSuite(PsfFluxTransformTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
