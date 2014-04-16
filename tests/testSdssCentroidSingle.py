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

numpy.random.seed(1234)

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each measRecord within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")

class SFMTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        catalog, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox)
        MakeTestData.fillImages(catalog, exposure)
        catalog.writeFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "calexp-0A.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "ref-0A.fits"))
    

    def tearDown(self):
        os.unlink(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "calexp-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "ref-0A.fits"))

    def testAlgorithm(self):

        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        exposure = lsst.afw.image.ExposureF(path)
        #  catalog with footprints, but not measurement fields added
        path = os.path.join(DATA_DIR, 'truthcat-0A.fits')
        srccat = SourceCatalog.readFits(path)
        #  catalog with footprints, but no measurement fields added
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in srccat}
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        replaced = lsst.afw.image.ExposureF(path)
        noiseReplacer = NoiseReplacer(replaced, footprints, sfm_config.noiseSource,
                          sfm_config.noiseOffset, sfm_config.noiseSeed)
        
        # add the measurement fields to the outputSchema and make a catalog with it
        # then extend with the mapper to copy the extant data
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        outschema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        sfm_config.plugins = ["centroid.peak", "base_SdssCentroid"]
        sfm_config.slots.centroid = "centroid.peak"
        sfm_config.slots.shape = None
        sfm_config.slots.psfFlux = None
        sfm_config.slots.modelFlux = None
        sfm_config.slots.apFlux = None
        sfm_config.slots.instFlux = None
        task = SingleFrameMeasurementTask(outschema, flags, config=sfm_config)
        measCat = SourceCatalog(outschema)
        measCat.extend(srccat, mapper=mapper)
        # now run the SFM task with the test plugin
        task.run(measCat, exposure)

        # The test plugin adds the footprint flux and the background (surrounding) flux
        # to the schema.  This test then loops through the sources and tries to produce
        # the same results
        mi = exposure.getMaskedImage()
        truthFluxkey = srccat.getSchema().find("truth.flux").key
        schema = measCat.getSchema()
        for i in range(len(measCat)):
            record = measCat[i]
            x = record.get("base_SdssCentroid_x")
            y = record.get("base_SdssCentroid_y")
            centroid = record.getCentroid() 
            print x,y, centroid
            self.assertFalse(record.get("base_SdssCentroid_flag"))
            self.assertFalse(record.get("base_SdssCentroid_flag_badData"))
            self.assertFalse(record.get("base_SdssCentroid_flag_edge"))
            centroid = record.getCentroid() 
            self.assertClose(x, centroid.getX(), atol=None, rtol=.01)
            self.assertClose(y, centroid.getY(), atol=None, rtol=.01)
    


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
