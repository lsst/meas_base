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

class TestCentroidConfig(SingleFramePluginConfig):
    fractional = lsst.pex.config.Field(dtype=bool, default=True,
                    doc="whether to center on fractional pixels")


class TestCentroid(SingleFramePlugin):
    ConfigClass = TestCentroidConfig

    def __init__(self, config, name, schema=None, flags=None, others=None, metadata=None):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.xKey = schema.addField("centroid.x", type=float, doc="x component relative to image",
                                    units="pixels")
        self.yKey = schema.addField("centroid.y", type=float, doc="y component relative to image",
                                    units="pixels")
        self.cKey = schema.find("centroid.sdss").key

    def measure(self, measRecord, exposure):
        c = measRecord.get(self.ckey)
        measRecord.set(self.xKey, c.getX())
        measRecord.set(self.yKey, c.getY())

class TestFluxConfig(SingleFramePluginConfig):
    pass

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each measRecord within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class TestFlux(SingleFramePlugin):
    ConfigClass = TestFluxConfig

    def __init__(self, config, name, schema=None, flags=None, others=None, metadata=None):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.fluxKey = schema.addField("test.flux", type=float, doc="sum of flux in object footprint",
                                       units="dn")
        self.fluxCountKey = schema.addField("test.fluxCount", type=int,
                                            doc="number of pixels in object footprint", units="pixels^2")
        self.backKey = schema.addField("test.back", type=float, doc="avg of flux in background", units="dn")
        self.backCountKey = schema.addField("test.backCount", type=int,
                                            doc="number of pixels in surrounding background",
                                            units="pixels^2")

    def measure(self, measRecord, exposure):
        foot = measRecord.getFootprint()
        image = exposure.getMaskedImage().getImage()
        array = image.getArray()

        # sum the footprint area for this measRecord
        sumArray = numpy.ndarray((foot.getArea()), array.dtype)
        lsst.afw.detection.flattenArray(foot, image.getArray(), sumArray, image.getXY0())
        flux = sumArray.sum(dtype=numpy.float64)
        area = foot.getArea()
        measRecord.set(self.fluxKey, flux)
        measRecord.set(self.fluxCountKey, area)

        # Now find an area which is 100 pixels larger in all directions than the foot.getBBox()
        fBBox = foot.getBBox()
        border = 100
        xmin = fBBox.getMinX() - border
        ymin = fBBox.getMinY() - border
        xmax = fBBox.getMaxX() + border
        ymax = fBBox.getMaxY() + border
        x0 = image.getX0()
        y0 = image.getY0()
        if xmin < x0: xmin = x0
        if ymin < y0: ymin = y0
        if xmax > (x0 + exposure.getWidth()): xmax = x0+exposure.getWidth()
        if ymax > (y0 + exposure.getHeight()): ymax = y0+exposure.getHeight()
        bigArraySub = array[ymin-y0:ymax-y0, xmin-x0:xmax-x0]
        bigFlux = bigArraySub.sum(dtype=numpy.float64)
        bigArea = (ymax-ymin)*(xmax-xmin)
        measRecord.set(self.backKey, bigFlux - flux)
        measRecord.set(self.backCountKey, bigArea - area)

SingleFramePlugin.registry.register("test.flux", TestFlux)

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")

class SFMTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        catalog, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox)
        MakeTestData.fillImages(catalog, exposure)
        catalog.writeFits(os.path.join(DATA_DIR, "truthcat-0C.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "calexp-0C.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "ref-0C.fits"))
    

    def tearDown(self):
        os.unlink(os.path.join(DATA_DIR, "truthcat-0C.fits"))
        os.unlink(os.path.join(DATA_DIR, "calexp-0C.fits"))
        os.unlink(os.path.join(DATA_DIR, "ref-0C.fits"))

    #  Run the measurement (sfm) task with its default plugins.  Any run to completion is successful
    def testRunMeasurement(self):
        path = os.path.join(DATA_DIR, 'calexp-0C.fits')
        exposure = lsst.afw.image.ExposureF(path)
        #  catalog with footprints, but not measurement fields added
        path = os.path.join(DATA_DIR, 'truthcat-0C.fits')
        srccat = SourceCatalog.readFits(path)
        flags = MeasurementDataFlags()

        # Read a catalog which should be the same as the catalog of processCcd
        # prior to measurement.  Create an empty catalog with the same schema
        # plus the schema items for the SFM task, then transfer the existing data
        # to the new catalog
        schema = srccat.getSchema()
        config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        task = SingleFrameMeasurementTask(schema, flags, config=config)
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        schema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        config.plugins = ["centroid.peak"]
        config.slots.centroid = None #"centroid.peak"
        config.slots.shape = None
        config.slots.psfFlux = None 
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.instFlux = None
        task = SingleFrameMeasurementTask(schema, flags, config=config)
        measCat = SourceCatalog(schema)
        measCat.extend(srccat, mapper=mapper)

        # Then run the default SFM task.  Results not checked
        task.run(measCat, exposure)

    #  This test really tests both that a plugin can measure things correctly,
    #  and that the noise replacement mechanism works in situ.
    #  The test uses the same replacement image as the test above, with the
    #  default NoiseReplacer seed.  This test just checks to be sure that the
    #  base.py replacement mechanism is still working
    def testFluxPlugin(self):

        print "testFluxPlugin"
        path = os.path.join(DATA_DIR, 'calexp-0C.fits')
        exposure = lsst.afw.image.ExposureF(path)
        #  catalog with footprints, but not measurement fields added
        path = os.path.join(DATA_DIR, 'truthcat-0C.fits')
        srccat = SourceCatalog.readFits(path)
        #  catalog with footprints, but no measurement fields added
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in srccat}
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        path = os.path.join(DATA_DIR, 'calexp-0C.fits')
        replaced = lsst.afw.image.ExposureF(path)
        noiseReplacer = NoiseReplacer(sfm_config.noiseReplacer, replaced, footprints)

        # add the measurement fields to the outputSchema and make a catalog with it
        # then extend with the mapper to copy the extant data
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        outschema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        sfm_config.plugins = ["centroid.peak", "test.flux"]
        sfm_config.slots.centroid = None #"centroid.peak"
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
        schema = measCat.getSchema()
        fluxkey = schema.find("test.flux").key
        fluxCountKey = schema.find("test.fluxCount").key
        backKey = schema.find("test.back").key
        backCountKey = schema.find("test.backCount").key
        truthFluxkey = srccat.getSchema().find("truth_flux").key
   
        # Test all the records to be sure that the measurement mechanism works for total flux
        # And that the area surrounding the footprint has the expected replacement pixels
        for i in range(len(measCat)):
            record = measCat[i]
            # First check to be sure that the flux measured by the plug-in is correct
            # This repeats the algorithm used by the NoiseReplacer, so not sure how useful it is
            foot = record.getFootprint()
            if foot.isHeavy():
                heavy = lsst.afw.detection.HeavyFootprintF.cast(foot)
            else:
                heavy = lsst.afw.detection.makeHeavyFootprint(foot, mi)
            sum = heavy.getImageArray().sum()
            count = foot.getArea()
            # get the values produced by the plugin
            flux = record.get(fluxkey)
            fluxCount = record.get(fluxCountKey)
            # Test 1:  the flux in the footprint area should measure the same as during the SFM run
            self.assertEqual(count,fluxCount)
            truthFlux = srccat[i].get(truthFluxkey)

            #  The flux reported by the test.flux plugin should be VERY close to the footprint sum
            #      the two differing only because of noise
            #  The flux reported by the test.flux plugin should be close to the truthFlux, but could
            #      differ due to finite aperature effects. 
            self.assertClose(sum, flux, atol=None, rtol=.0001)
            self.assertClose(sum, truthFlux, atol=None, rtol=.03)


            # Now find an area which is 100 pixels larger in all directions than the foot.getBBox()
            # Measure it against the replaced.fits image to be sure no extra flux is appearing
            # due to incorrect background fill in of any of the other objects
            #  heavy.insert(replaced)
            fbbox = foot.getBBox()
            border = 100

            x0 = mi.getImage().getX0()
            y0 = mi.getImage().getY0()
            xmin = fbbox.getMinX() - border
            ymin = fbbox.getMinY() - border
            xmax = fbbox.getMaxX() + border
            ymax = fbbox.getMaxY() + border
            if xmin < x0: xmin = x0
            if ymin < y0: ymin = y0
            if xmax > (exposure.getWidth()+x0): xmax = exposure.getWidth()+y0
            if ymax > (exposure.getHeight()+y0): ymax = exposure.getHeight()+y0

            # get the sum of the entire bordered area
            arraysub = replaced.getMaskedImage().getImage().getArray()[ymin-y0:ymax-y0, xmin-x0:xmax-x0]
            bigFlux = arraysub.sum(dtype=numpy.float64)

            # get the sum of the footprint area
            sumarray = numpy.ndarray((foot.getArea()), replaced.getMaskedImage().getImage().getArray().dtype)
            lsst.afw.detection.flattenArray(foot, replaced.getMaskedImage().getImage().getArray(),
                sumarray, replaced.getMaskedImage().getXY0())
            repFlux = sumarray.sum(dtype=numpy.float64)

            newBackCount = (ymax-ymin)*(xmax-xmin) - count
            newBack = bigFlux - repFlux
            back = record.get(backKey)
            backCount = record.get(backCountKey)


            # Test 2:  the area surrounding the object should be the same as what was measured
            #          during the actual run
            self.assertEqual(backCount, newBackCount)
            #self.assertClose(back,newBack, atol=None, rtol=.03)


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
