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
import lsst.meas.base

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data")

# n.b. Some tests here depend on the noise realization in the test data;
# for the one currently checked, they passes, but it may not if the test
# data is regenerated.

class PsfFluxTestCase(lsst.utils.tests.TestCase):
    """Test case for the PsfFlux algorithm/plugins
    """

    def setUp(self):
        self.truth = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-01.fits"))
        self.calexp = lsst.afw.image.ExposureF.readFits(os.path.join(DATA_DIR, "calexp-01.fits"))
        self.record = self.truth[0]
        self.fluxKey = self.truth.schema.find("truth.flux").key
        self.centroidKey = self.truth.schema.find("truth.centroid").key

    def tearDown(self):
        del self.truth
        del self.calexp
        del self.record
        del self.fluxKey
        del self.centroidKey

    def testOverloads(self):
        """Verify that we get reasonable answers when we call apply() directly on an isolated source,
        and that these results are identical regardless of which overload we use.
        """
        result1 = lsst.meas.base.PsfFluxAlgorithm.apply(
            self.calexp, self.record.getFootprint(),
            self.record.get(self.centroidKey)
            )
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
        self.assertClose(result1.flux, self.record.get(self.fluxKey), atol=result1.fluxSigma)

    def testMasking(self):
        badPoint = lsst.afw.geom.Point2I(self.record.get(self.centroidKey)) + lsst.afw.geom.Extent2I(3,4)
        imageArray = self.calexp.getMaskedImage().getImage().getArray()
        maskArray = self.calexp.getMaskedImage().getMask().getArray()
        badMask = self.calexp.getMaskedImage().getMask().getPlaneBitMask("BAD")
        imageArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] = numpy.inf
        maskArray[badPoint.getY() - self.calexp.getY0(), badPoint.getX() - self.calexp.getX0()] |= badMask
        # Should get an infinite value exception, because we didn't mask that one pixel
        lsst.utils.tests.assertRaisesLsstCpp(
            self, lsst.meas.base.PixelValueError,
            lsst.meas.base.PsfFluxAlgorithm.apply,
            self.calexp, self.record.getFootprint(),
            self.record.get(self.centroidKey)
            )
        # If we do mask it, we should get a reasonable result
        ctrl = lsst.meas.base.PsfFluxAlgorithm.Control()
        ctrl.badMaskPlanes = ["BAD"]
        result = lsst.meas.base.PsfFluxAlgorithm.apply(
            self.calexp, self.record.getFootprint(),
            self.record.get(self.centroidKey),
            ctrl
            )
        self.assertClose(result.flux, self.record.get(self.fluxKey), atol=result.fluxSigma)
        # If we mask the whole image, we should get a MeasurementError
        maskArray[:,:] |= badMask
        lsst.utils.tests.assertRaisesLsstCpp(
            self, lsst.meas.base.MeasurementError,
            lsst.meas.base.PsfFluxAlgorithm.apply,
            self.calexp, self.record.getFootprint(),
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
            subExposure, self.record.getFootprint(),
            self.record.get(self.centroidKey),
            )
        self.assertClose(result.flux, self.record.get(self.fluxKey), atol=result.fluxSigma)
        self.assertTrue(result.getFlag(lsst.meas.base.PsfFluxAlgorithm.EDGE))

    def testSubImage(self):
        """Test that we raise MeasurementError when there's no PSF.
        """
        self.calexp.setPsf(None)
        lsst.utils.tests.assertRaisesLsstCpp(
            self, lsst.meas.base.MeasurementError,
            lsst.meas.base.PsfFluxAlgorithm.apply,
            self.calexp, self.record.getFootprint(),
            self.record.get(self.centroidKey),
            )

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
