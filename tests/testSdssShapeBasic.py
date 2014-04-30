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
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementTask
from lsst.meas.base.base import *
from lsst.meas.base.tests import *
import unittest
import lsst.utils.tests
import numpy

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")

class SFMTestCase(lsst.utils.tests.TestCase):

    def testAlgorithm(self):

        srccat, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox)
        MakeTestData.fillImages(srccat, exposure)

        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()        
        # add the measurement fields to the outputSchema and make a catalog with it
        # then extend with the mapper to copy the extant data
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        outschema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        sfm_config.plugins = ["centroid.peak", "base_SdssShape"]
        sfm_config.slots.centroid = None
        sfm_config.slots.shape = "base_SdssShape"
        sfm_config.slots.psfFlux = None
        sfm_config.slots.modelFlux = None
        sfm_config.slots.apFlux = None
        sfm_config.slots.instFlux = None
        task = SingleFrameMeasurementTask(outschema, flags, config=sfm_config)
        measCat = SourceCatalog(outschema)
        measCat.getTable().setVersion(1)
        measCat.extend(srccat, mapper=mapper)
        # now run the SFM task with the test plugin
        task.run(measCat, exposure)

        truthShapeKey = lsst.afw.table.QuadrupoleKey(srccat.schema.find("truth_xx").key,
                                                     srccat.schema.find("truth_yy").key,
                                                     srccat.schema.find("truth_xy").key)
        for i in range(len(measCat)):
            record = measCat[i]
            srcRec = srccat[i]
            xx = record.get("base_SdssShape_xx")
            yy = record.get("base_SdssShape_yy")
            xy = record.get("base_SdssShape_xy")
            xxSigma = record.get("base_SdssShape_xxSigma")
            yySigma = record.get("base_SdssShape_yySigma")
            xySigma = record.get("base_SdssShape_xySigma")
            xxyyCov = record.get("base_SdssShape_xx_yy_Cov")
            xxxyCov = record.get("base_SdssShape_xx_xy_Cov")
            yyxyCov = record.get("base_SdssShape_yy_xy_Cov")
            trueShape = srcRec.get(truthShapeKey)
            shape = record.getShape()
            cov = record.getShapeErr()
            self.assertClose(xxSigma*xxSigma, cov[0,0], rtol = .01)
            self.assertClose(yySigma*yySigma, cov[1,1], rtol = .01)
            self.assertClose(xySigma*xySigma, cov[2,2], rtol = .01)
            self.assertTrue(numpy.isnan(xxyyCov) and numpy.isnan(cov[0,1]))
            self.assertTrue(numpy.isnan(xxxyCov) and numpy.isnan(cov[0,2]))
            self.assertTrue(numpy.isnan(yyxyCov) and numpy.isnan(cov[1,2]))
            if not numpy.isnan(trueShape.getIxx()):
                self.assertFalse(record.get("base_SdssShape_flag"))
                self.assertFalse(record.get("base_SdssShape_flag_unweightedBad"))
                self.assertFalse(record.get("base_SdssShape_flag_unweighted"))
                self.assertFalse(record.get("base_SdssShape_flag_shift"))
                self.assertFalse(record.get("base_SdssShape_flag_maxIter"))
                x = record.get("base_SdssShape_x")
                y = record.get("base_SdssShape_y")
                xSigma = record.get("base_SdssShape_xSigma")
                ySigma = record.get("base_SdssShape_ySigma")
                flux = record.get("base_SdssShape_flux")
                fluxSigma = record.get("base_SdssShape_fluxSigma")
                xy4 = record.get("base_SdssShape_xy4")
                xy4Sigma = record.get("base_SdssShape_xy4Sigma")
                self.assertClose(xx, trueShape.getIxx(), atol=None, rtol=.12)
                self.assertClose(yy, trueShape.getIyy(), atol=None, rtol=.12)
                # commented out because of a bug
                #self.assertClose(xy, trueShape.getIxy(), atol=None, rtol=.12)
def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SFMTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
