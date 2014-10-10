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

import math
import os
from lsst.afw.table import Schema, SchemaMapper, SourceCatalog, SourceTable
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementTask
from lsst.meas.base.base import *
from lsst.meas.base.tests import *
import unittest
import lsst.utils.tests
import numpy

numpy.random.seed(1234)

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")

class SdssCentroidTestCase(lsst.meas.base.tests.AlgorithmTestCase):

    def setUp(self):
        lsst.meas.base.tests.AlgorithmTestCase.setUp(self)
        sfmConfig = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        # add the measurement fields to the outputSchema and make a catalog with it
        # then extend with the mapper to copy the extant data
        mapper = SchemaMapper(self.truth.getSchema())
        mapper.addMinimalSchema(self.truth.getSchema())
        outSchema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        sfmConfig.plugins = ["base_PeakCentroid", "base_SdssCentroid"]
        sfmConfig.slots.centroid = "base_SdssCentroid"
        sfmConfig.slots.shape = None
        sfmConfig.slots.psfFlux = None
        sfmConfig.slots.modelFlux = None
        sfmConfig.slots.apFlux = None
        sfmConfig.slots.instFlux = None
        self.task = SingleFrameMeasurementTask(outSchema, flags, config=sfmConfig)
        self.measCat = SourceCatalog(outSchema)
        self.measCat.defineCentroid("base_SdssCentroid")
        self.measCat.extend(self.truth, mapper=mapper)

    def tearDown(self):
        lsst.meas.base.tests.AlgorithmTestCase.tearDown(self)
        del self.task
        del self.measCat

    def testAlgorithm(self):
        # now run the SFM task with the test plugin
        self.task.run(self.measCat, self.calexp)

        truthFluxkey = self.truth.getSchema().find("truth_flux").key
        for i in range(len(self.measCat)):
            record = self.measCat[i]
            centroid = record.getCentroid()
            cov = record.getCentroidErr()
            peakX = record.get("base_PeakCentroid_x")
            peakY = record.get("base_PeakCentroid_y")
            x = record.get("base_SdssCentroid_x")
            y = record.get("base_SdssCentroid_y")
            xerr = record.get("base_SdssCentroid_xSigma")
            yerr = record.get("base_SdssCentroid_ySigma")
            xycov = record.get("base_SdssCentroid_x_y_Cov")

            self.assertFalse(record.get("base_SdssCentroid_flag"))
            self.assertFalse(record.get("base_SdssCentroid_flag_badData"))
            self.assertFalse(record.get("base_SdssCentroid_flag_edge"))
            self.assertClose(peakX, x, atol=None, rtol=.02)
            self.assertClose(peakY, y, atol=None, rtol=.02)

            self.assertEqual(x, record.getCentroid().getX())
            self.assertEqual(y, record.getCentroid().getY())
            self.assertClose(xerr*xerr, cov[0,0], rtol=.01)
            self.assertClose(yerr*yerr, cov[1,1], rtol=.01)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SdssCentroidTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
