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

import lsst.utils.tests
import lsst.afw.detection
import lsst.afw.table
import lsst.meas.base.tests

numpy.random.seed(1234)

@lsst.meas.base.register("test_NoiseReplacer")
class NoiseReplacerTestPlugin(lsst.meas.base.SingleFramePlugin):
    """A measurement plugin that simply sums flux inside and outside the source's footprint."""

    @staticmethod
    def getExecutionOrder():
        return 2.0

    def __init__(self, config, name, schema, metadata):
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.insideKey = schema.addField("%s_inside" % name, type=float, doc="flux inside footprint")
        self.outsideKey = schema.addField("%s_outside" % name, type=float, doc="flux outside footprint")

    def measure(self, measRecord, exposure):
        footprint = measRecord.getFootprint()
        fullArray = exposure.getMaskedImage().getImage().getArray()
        insideArray = numpy.zeros(footprint.getArea(), dtype=fullArray.dtype)
        lsst.afw.detection.flattenArray(footprint, fullArray, insideArray, exposure.getXY0())
        insideFlux = float(insideArray.sum())
        outsideFlux = float(fullArray.sum()) - insideFlux
        measRecord.set(self.insideKey, insideFlux)
        measRecord.set(self.outsideKey, outsideFlux)


class NoiseReplacerTestCase(lsst.meas.base.tests.AlgorithmTestCase):

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(240, 260))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(120000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))
        with self.dataset.addBlend() as family:
            family.addChild(110000.0, lsst.afw.geom.Point2D(65.2, 150.7),
                            lsst.afw.geom.ellipses.Quadrupole(7, 5, -1))
            family.addChild(140000.0, lsst.afw.geom.Point2D(72.3, 149.1))
            family.addChild(90000.0, lsst.afw.geom.Point2D(68.5, 156.9))

    def testSingleFrameMeasurement(self):
        """Test that replacing sources with noise works as used in SingleFrameMeasurementTask,
        by comparing flux inside and outside source Footprints on an extremely high S/N image."""
        task = self.makeSingleFrameMeasurementTask("test_NoiseReplacer")
        exposure, catalog = self.dataset.realize(1.0, task.schema)
        task.run(exposure, catalog)
        sumVariance = exposure.getMaskedImage().getVariance().getArray().sum()
        for record in catalog:
            self.assertClose(record.get("test_NoiseReplacer_inside"), record.get("truth_flux"), rtol=1E-3)
            # n.b. Next line checks that a random value is correct to a statistical 1-sigma prediction;
            # some RNG seeds may cause it to fail (indeed, 67% should)
            self.assertLess(record.get("test_NoiseReplacer_outside"), numpy.sqrt(sumVariance))

    def tearDown(self):
        del self.bbox
        del self.dataset

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(NoiseReplacerTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
