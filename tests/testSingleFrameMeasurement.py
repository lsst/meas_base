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
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import *
from lsst.meas.base.base import *
import unittest
import lsst.utils.tests
import numpy

class TestCentroidConfig(SingleFramePluginConfig):
    fractional = lsst.pex.config.Field(dtype=bool, default=True,
                    doc="whether to center on fractional pixels")


class TestCentroid(SingleFramePlugin):

    ConfigClass = TestCentroidConfig
    doMeasureSingle = True
    doMeasureMulti = False
    def __init__(self, config, name, schema=None, flags=None, others=None, metadata=None):

            schema.addField("centroid.x", type=float, doc="x component relative to image", units="pixels")
            schema.addField("centroid.y", type=float, doc="x component relative to image", units="pixels")

    def measureSingle(self, exposure, source):
        schema = source.getSchema()
        ckey = schema.find("centroid.sdss").key
        ykey = schema.find("centroid.y").key
        xkey = schema.find("centroid.x").key
        c = source.get(ckey)
        source.set(xkey, c.getX())
        source.set(ykey, c.getY())
        return

    def measureMulti(self, exposure, sources):
        return

class TestFluxConfig(SingleFramePluginConfig):
    pass

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each source within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class TestFlux(SingleFramePlugin):
    ConfigClass = TestFluxConfig
    doMeasureSingle = True
    doMeasureMulti = False

    def __init__(self, config, name, schema=None, flags=None, others=None, metadata=None):

        schema.addField("test.flux", type=float, doc="sum of flux in object footprint", units="")
        schema.addField("test.fluxcount", type=int, doc="number of pixels in object footprint", units="")
        schema.addField("test.back", type=float, doc="avg of flux in background", units="")
        schema.addField("test.backcount", type=int, doc="number of pixels in surrounding background",
            units="")
        self.config = config

    def measureSingle(self, exposure, source):

        schema = source.getSchema()
        fluxkey = schema.find("test.flux").key
        fluxcountkey = schema.find("test.fluxcount").key
        backkey = schema.find("test.back").key
        backcountkey = schema.find("test.backcount").key
        id = source.getId()
        foot = source.getFootprint()
        image = exposure.getMaskedImage().getImage()
        array = image.getArray()

        # sum the footprint area for this source
        sumarray = numpy.ndarray((foot.getArea()), array.dtype)
        lsst.afw.detection.flattenArray(foot, image.getArray(), sumarray, image.getXY0())
        flux = sumarray.sum(dtype=numpy.float64)
        area = foot.getArea()
        source.set(fluxkey, flux)
        source.set(fluxcountkey, area)

        # Now find an area which is 100 pixels larger in all directions than the foot.getBBox()
        fbbox = foot.getBBox()
        border = 100
        xmin = fbbox.getMinX() - border
        ymin = fbbox.getMinY() - border
        xmax = fbbox.getMaxX() + border
        ymax = fbbox.getMaxY() + border
        x0 = image.getX0()
        y0 = image.getY0()
        if xmin < x0: xmin = x0
        if ymin < y0: ymin = y0
        if xmax > (x0 + exposure.getWidth()): xmax = x0+exposure.getWidth()
        if ymax > (y0 + exposure.getHeight()): ymax = y0+exposure.getHeight()
        bigarraysub = array[ymin-y0:ymax-y0, xmin-x0:xmax-x0]
        bigflux = bigarraysub.sum(dtype=numpy.float64)
        bigarea = (ymax-ymin)*(xmax-xmin)
        source.set(backkey, bigflux - flux)
        source.set(backcountkey, bigarea - area)


    def measureMulti(self, exposure, sources):
        return

SingleFramePlugin.registry.register("test.flux", TestFlux)

class SFMTestCase(lsst.utils.tests.TestCase):
    # Test the Noise Replacement mechanism.  This is an extremely cursory test, just
    # to be sure that the patterns used to fill in all of the footprints are more or
    # less random and Gaussian.  Probably should have a really normality test here
    def testANoiseReplacement(self):
        exposure = lsst.afw.image.ExposureF("data/exposure.fits.gz")
        #  catalog with footprints, but not measurement fields added
        srccat = SourceCatalog.readFits("data/measCatNull.fits.gz")
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in srccat}
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
     
        # create an exposure which is identical to the exposure created in actual measurement runs
        # this has random noise in place of the source footprints
        replaced = lsst.afw.image.ExposureF("data/exposure.fits.gz")
        noiseReplacer = NoiseReplacer(replaced, footprints, sfm_config.noiseSource,
                          sfm_config.noiseOffset, sfm_config.noiseSeed)

        # accuulate the variation of the mean for each filled region in std units
        normtest = numpy.ndarray(len(srccat), numpy.float32)
        for i in range(len(srccat)):
            record = srccat[i]
            # First check to be sure that the flux measured by the plug-in is correct
            # This repeats the algorithm used by the NoiseReplacer, so not sure how useful it is
            foot = record.getFootprint()
            # get the sum of the footprint area
            sumarray = numpy.ndarray((foot.getArea()), replaced.getMaskedImage().getImage().getArray().dtype)
            lsst.afw.detection.flattenArray(foot, replaced.getMaskedImage().getImage().getArray(),
                sumarray, replaced.getMaskedImage().getXY0())
            repflux = sumarray.sum(dtype=numpy.float64)

            # Test: the fill in of the area should be consistent with the noiseReplacer
            # mean=0.702828, std=27.4483
            noisemean = noiseReplacer.noiseGenMean
            noisestd = noiseReplacer.noiseGenStd
            if noisemean == None or noisestd == None:
                return   # if the value is not available, skip this test
            popstd = math.sqrt((noisestd*noisestd)/len(sumarray))
            normtest[i] = (sumarray.mean()-noisemean)/popstd

        #  These values pass using the current noise generator.  Since the seed is fixed, it will always pass.
        #  Probably should put a random seed and normality test in place of this.
        self.assertTrue(normtest.min() > -3.5)
        self.assertTrue(normtest.max() < 3.5)
        self.assertTrue(abs(normtest.mean()) < .02)
        self.assertTrue(abs(normtest.std() - 1.0) < .02)

    #  Run the measurement (sfm) task with its default plugins.  Any run to completion is successful
    def testRunMeasurement(self):

        print "testRunMeasurement"
        exposure = lsst.afw.image.ExposureF("data/exposure.fits.gz")
        flags = MeasurementDataFlags()

        # Read a catalog which should be the same as the catalog of processCcd
        # prior to measurement.  Create an empty catalog with the same schema
        # plus the schema items for the SFM task, then transfer the existing data
        # to the new catalog
        srccat = SourceCatalog.readFits("data/measCatNull.fits.gz")
        schema = srccat.getSchema()
        config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        task = SingleFrameMeasurementTask(schema, flags, config=config)
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        schema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        task = SingleFrameMeasurementTask(schema, flags)
        measCat = SourceCatalog(schema)
        measCat.extend(srccat, mapper=mapper)

        # Then run the default SFM task.  Results not checked
        task.run(exposure, measCat)
        measCat.writeFits("test.fits")

    #  This test really tests both that a plugin can measure things correctly,
    #  and that the noise replacement mechanism works in situ.
    #  The test uses the same replacement image as the test above, with the
    #  default NoiseReplacer seed.  This test just checks to be sure that the
    #  base.py replacement mechanism is still working

    def testFluxPlugin(self):

        exposure = lsst.afw.image.ExposureF("data/exposure.fits.gz")
        #  catalog with footprints, but no measurement fields added
        srccat = SourceCatalog.readFits("data/measCatNull.fits.gz")
        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in srccat}
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        sfm_config.plugins.names.add("test.flux")
        replaced = lsst.afw.image.ExposureF("data/exposure.fits.gz")
        noiseReplacer = NoiseReplacer(replaced, footprints, sfm_config.noiseSource,
                          sfm_config.noiseOffset, sfm_config.noiseSeed)

        # add the measurement fields to the outputSchema and mack a catalog with it
        # then extend with the mapper to copy the extant data
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        outschema = mapper.getOutputSchema()
        flags = MeasurementDataFlags()
        task = SingleFrameMeasurementTask(outschema, flags, config=sfm_config)
        measCat = SourceCatalog(outschema)
        measCat.extend(srccat, mapper=mapper)

        # now run the SFM task with the test plugin
        task.run(exposure, measCat)

        # The test plugin adds the footprint flux and the background (surrounding) flux
        # to the schema.  This test then loops through the sources and tries to produce
        # the same results
        mi = exposure.getMaskedImage()
        schema = measCat.getSchema()
        fluxkey = schema.find("test.flux").key
        fluxcountkey = schema.find("test.fluxcount").key
        backkey = schema.find("test.back").key
        backcountkey = schema.find("test.backcount").key

        # Test all the records to be sure that the measurement mechanism works for total flux
        # And that the area surrounding the footprint has the expected replacement pixels
        for i in range(len(measCat)):
            record = measCat[i]
            # First check to be sure that the flux measured by the plug-in is correct
            # This repeats the algorithm used by the NoiseReplacer, so not sure how useful it is
            foot = record.getFootprint()
            if foot.isHeavy():
                heavy = afwDet.cast_HeavyFootprintF(foot)
            else:
                heavy = lsst.afw.detection.makeHeavyFootprint(foot, mi)
            noise = lsst.afw.detection.makeHeavyFootprint(foot, replaced.getMaskedImage())
            sumarray = numpy.ndarray((foot.getArea()), mi.getImage().getArray().dtype)
            lsst.afw.detection.flattenArray(foot, mi.getImage().getArray(), sumarray, mi.getImage().getXY0())
            sum = sumarray.sum(dtype=numpy.float64)
            count = foot.getArea()
            # get the values produced by the plugin
            flux = record.get(fluxkey)
            fluxcount = record.get(fluxcountkey)
            # Test 1:  the flux in the footprint area should measure the same as during the SFM run
            self.assertEqual(count,fluxcount)
            self.assertEqual(sum,flux)


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
            bigflux = arraysub.sum(dtype=numpy.float64)

            # get the sum of the footprint area
            sumarray = numpy.ndarray((foot.getArea()), replaced.getMaskedImage().getImage().getArray().dtype)
            lsst.afw.detection.flattenArray(foot, replaced.getMaskedImage().getImage().getArray(),
                sumarray, replaced.getMaskedImage().getXY0())
            repflux = sumarray.sum(dtype=numpy.float64)

            newbackcount = (ymax-ymin)*(xmax-xmin) - count
            newback = bigflux - repflux
            back = record.get(backkey)
            backcount = record.get(backcountkey)


            # Test 2:  the area surrounding the object should be the same as what was measured
            #          during the actual run
            self.assertEqual(backcount, newbackcount)
            # dont know self.assertEqual(back,newback)


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
