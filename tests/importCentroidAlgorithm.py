#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import re
import os
import glob
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDetection
import lsst.meas.base
import lsst.meas.algorithms
import lsst.utils.tests as utilsTests

import testlib

try:
    type(verbose)
except NameError:
    verbose = 0

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CentroidTestCase(unittest.TestCase):
    """A test case for centroiding"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testApplyCentroid(self):
        """Test that we can instantiate and play with SillyMeasureCentroid"""

        for imageFactory in (
                             afwImage.MaskedImageF,
                             afwImage.MaskedImageD,
                             ):
            im = imageFactory(afwGeom.ExtentI(100, 100))
            exp = afwImage.makeExposure(im)
            for offset in (0,1,2):
                control = testlib.SillyCentroidControl()
                control.param = offset
                results = testlib.SillyCentroidAlgorithm.Result()
                x, y = 10, 20
                testlib.SillyCentroidAlgorithm.apply( exp, afwGeom.Point2D(x, y), results, control)
                self.assertEqual(x, results.getCentroid().getX() - offset)
                self.assertEqual(y, results.getCentroid().getY() - offset)


    def testMeasureCentroid(self):
        """Test that we can use our silly centroid through the usual Tasks"""
        lsst.meas.base.WrappedSingleFramePlugin.generate(testlib.SillyCentroidAlgorithm, executionOrder=0.0)
        x, y = 10, 20

        im = afwImage.MaskedImageF(afwGeom.ExtentI(512, 512))
        im.set(0)
        arr = im.getImage().getArray()
        arr[y,x] = 1
        exp = afwImage.makeExposure(im)

        detConfig = lsst.meas.algorithms.SourceDetectionConfig()
        detConfig.thresholdValue = 0.5
        detConfig.thresholdType = "value"
        schema = afwTable.SourceTable.makeMinimalSchema()
        det = lsst.meas.algorithms.SourceDetectionTask(schema=schema, config=detConfig)
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        flags = lsst.meas.base.MeasurementDataFlags()
        sfm_config.plugins = ["testlib_SillyCentroid"]
        sfm_config.plugins["testlib_SillyCentroid"].param = 5
        sfm_config.slots.centroid = "testlib_SillyCentroid"
        sfm_config.slots.shape = None
        sfm_config.slots.psfFlux = None
        sfm_config.slots.instFlux = None
        sfm_config.slots.apFlux = None
        sfm_config.slots.modelFlux = None
        task = lsst.meas.base.SingleFrameMeasurementTask(schema, flags, config=sfm_config)
        table = afwTable.SourceTable.make(schema)
        sources = det.makeSourceCatalog(table, exp, doSmooth=False, sigma=1.0).sources
        self.assertEqual(len(sources), 1)

        # now run the SFM task with the test plugin
        task.run(sources, exp)

        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].getY(), y + 5)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CentroidTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
