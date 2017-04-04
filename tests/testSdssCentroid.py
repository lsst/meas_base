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

from __future__ import absolute_import, division, print_function
from builtins import range
import unittest

import numpy as np


from lsst.meas.base.tests import (AlgorithmTestCase, CentroidTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)
import lsst.utils.tests

# n.b. Some tests here depend on the noise realization in the test data
# or from the numpy random number generator.
# For the current test data and seed value, they pass, but they may not
# if the test data is regenerated or the seed value changes.  I've marked
# these with an "rng dependent" comment.  In most cases, they test that
# the measured flux lies within 2 sigma of the correct value, which we
# should expect to fail sometimes.


class SdssCentroidTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.afw.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):

        del self.center
        del self.bbox
        del self.dataset

    def makeAlgorithm(self, ctrl=None):
        """Construct an algorithm (finishing a schema in the process), and return both.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.SdssCentroidControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.SdssCentroidAlgorithm(ctrl, "base_SdssCentroid", schema)
        return algorithm, schema

    def testSingleFramePlugin(self):
        """Test that we can call the algorithm through the SFM plugin interface."""
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(catalog, exposure)
        record = catalog[0]
        self.assertFalse(record.get("base_SdssCentroid_flag"))
        self.assertFalse(record.get("base_SdssCentroid_flag_edge"))
        self.assertClose(record.get("base_SdssCentroid_x"), record.get("truth_x"), rtol=0.005)
        self.assertClose(record.get("base_SdssCentroid_y"), record.get("truth_y"), rtol=0.005)

    def testMonteCarlo(self):
        """Test that we get exactly the right answer on an ideal sim with no noise, and that
        the reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        algorithm, schema = self.makeAlgorithm()
        exposure, catalog = self.dataset.realize(0.0, schema)
        record = catalog[0]
        x = record.get("truth_x")
        y = record.get("truth_y")
        flux = record.get("truth_flux")
        algorithm.measure(record, exposure)
        self.assertClose(record.get("base_SdssCentroid_x"), x, rtol=1E-4)
        self.assertClose(record.get("base_SdssCentroid_y"), y, rtol=1E-4)
        for noise in (0.001, 0.01):
            xList = []
            yList = []
            xSigmaList = []
            ySigmaList = []
            nSamples = 1000
            for repeat in range(nSamples):
                exposure, catalog = self.dataset.realize(noise*flux, schema)
                record = catalog[0]
                algorithm.measure(record, exposure)
                xList.append(record.get("base_SdssCentroid_x"))
                yList.append(record.get("base_SdssCentroid_y"))
                xSigmaList.append(record.get("base_SdssCentroid_xSigma"))
                ySigmaList.append(record.get("base_SdssCentroid_ySigma"))
            xMean = np.mean(xList)
            yMean = np.mean(yList)
            xSigmaMean = np.mean(xSigmaList)
            ySigmaMean = np.mean(ySigmaList)
            xStandardDeviation = np.std(xList)
            yStandardDeviation = np.std(yList)
            self.assertClose(xSigmaMean, xStandardDeviation, rtol=0.2)   # rng dependent
            self.assertClose(ySigmaMean, yStandardDeviation, rtol=0.2)   # rng dependent
            self.assertLess(xMean - x, 3.0*xSigmaMean / nSamples**0.5)   # rng dependent
            self.assertLess(yMean - y, 3.0*ySigmaMean / nSamples**0.5)   # rng dependent

    def testEdge(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        psfImage = exposure.getPsf().computeImage(self.center)
        # construct a box that won't fit the full PSF model
        bbox = psfImage.getBBox()
        bbox.grow(-5)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # we also need to install a smaller footprint, or NoiseReplacer complains before we even get to
        # measuring the centriod
        record = catalog[0]
        newFootprint = lsst.afw.detection.Footprint(bbox)
        peak = record.getFootprint().getPeaks()[0]
        newFootprint.addPeak(peak.getFx(), peak.getFy(), peak.getPeakValue())
        record.setFootprint(newFootprint)
        # just measure the one object we've prepared for
        task.measure(catalog, subImage)
        self.assertTrue(record.get("base_SdssCentroid_flag"))
        self.assertTrue(record.get("base_SdssCentroid_flag_edge"))

    def testNo2ndDerivative(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        # cutout a subimage around object in the test image
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(self.center), lsst.afw.geom.Extent2I(1, 1))
        bbox.grow(20)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # A completely flat image will trigger the no 2nd derivative error
        subImage.getMaskedImage().getImage().getArray()[:] = 0
        task.measure(catalog, subImage)
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag"))
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag_noSecondDerivative"))

    def testNotAtMaximum(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        # cutout a subimage around the object in the test image
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(self.center), lsst.afw.geom.Extent2I(1, 1))
        bbox.grow(20)
        subImage = lsst.afw.image.ExposureF(exposure, bbox)
        # zero out the central region, which will destroy the maximum
        subImage.getMaskedImage().getImage().getArray()[18:22, 18:22] = 0
        task.measure(catalog, subImage)
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag"))
        self.assertTrue(catalog[0].get("base_SdssCentroid_flag_notAtMaximum"))


class SdssCentroidTransformTestCase(CentroidTransformTestCase,
                                    SingleFramePluginTransformSetupHelper,
                                    lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.SdssCentroidControl
    algorithmClass = lsst.meas.base.SdssCentroidAlgorithm
    transformClass = lsst.meas.base.SdssCentroidTransform
    flagNames = ('flag', 'flag_edge', 'flag_badData')
    singleFramePlugins = ('base_SdssCentroid',)
    forcedPlugins = ('base_SdssCentroid',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
