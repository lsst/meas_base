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
import math
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import *
from lsst.meas.base.base import *
from lsst.daf.persistence.butler import *
import unittest
import lsst.utils.tests
import numpy

from lsst.meas.base.sfm import *
import lsst.afw.detection
import numpy

class TestCentroidConfig(SingleFrameAlgorithmConfig):
    fractional = lsst.pex.config.Field(dtype=bool, default=True,
                    doc="whether to center on fractional pixels")


class TestCentroid(SingleFrameAlgorithm):

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

class TestFluxConfig(SingleFrameAlgorithmConfig):
    pass

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each source within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class TestFlux(SingleFrameAlgorithm):
    ConfigClass = TestFluxConfig
    doMeasureSingle = True
    doMeasureMulti = False

    def __init__(self, config, name, schema=None, flags=None, others=None, metadata=None):

            schema.addField("test.flux", type=float, doc="sum of flux in object footprint", units="")
            schema.addField("test.fluxcount", type=int, doc="number of pixels in object footprint", units="")
            schema.addField("test.back", type=float, doc="avg of flux in background", units="")
            schema.addField("test.backcount", type=int, doc="number of pixels in surrounding background", units="")

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


SingleFrameAlgorithm.registry.register("test.flux", TestFlux)

class SFMTestCase(lsst.utils.tests.TestCase):

    #  Run the measurement (sfm) task with its default plugins
    #  Any run to completion is successful
    def testRunMeasurement(self):
        exposure = lsst.afw.image.ExposureF("data/exposure.fits")
        flags = MeasurementDataFlags()

        # Read a catalog which should be the same as the catalog of processCcd
        # prior to measurement.  Create an empty catalog with the same schema
        # plus the schema items for the SFM task, then transfer the existing data
        # to the new catalog
        srccat = SourceCatalog.readFits("data/measCatNull.fits")
        schema = srccat.getSchema()
        config = lsst.meas.base.sfm_algorithms.SingleFrameMeasurementConfig()
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
        
    #  This test really tests both that a plugin can measure things correctly,
    #  and that the noise replacement mechanism work.  The test assumes that
    #  the fits file replaced.fits is a correct noise replaced image with the
    #  default NoiseReplacer seed.  This test just checks to be sure that the
    #  base.py replacement mechanism is still working

    def testNoiseReplacement(self):
        exposure = lsst.afw.image.ExposureF("data/exposure.fits")
        replaced = lsst.afw.image.MaskedImageF("data/replaced.fits")
        config = lsst.meas.base.sfm_algorithms.SingleFrameMeasurementConfig()
        config.algorithms.names.add("test.flux")  
        flags = MeasurementDataFlags()

        #  catalog with footprints, but not measurement fields added
        srccat = SourceCatalog.readFits("data/measCatNull.fits")
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())

        # add the measurement fields to the outputSchema and mack a catalog with it
        # then extend with the mapper to copy the extant data
        outschema = mapper.getOutputSchema()
        task = SingleFrameMeasurementTask(outschema, flags, config=config)
        measCat = SourceCatalog(outschema)
        measCat.extend(srccat, mapper=mapper)

        # now run the SFM task with the test plugin
        task.run(exposure, measCat)

        # The test plugin adds the footprint flux and the background (surrounding) flux
        # to the schema.  This is not really a very interesting test, but we will do the
        # same thing here with the replaced.fits image and the heavy footprint for each
        # source to see if the SFM task gives the same results        
        mi = exposure.getMaskedImage()
        schema = measCat.getSchema()
        fluxkey = schema.find("test.flux").key
        fluxcountkey = schema.find("test.fluxcount").key
        backkey = schema.find("test.back").key
        backcountkey = schema.find("test.backcount").key
        normtest = numpy.ndarray(len(measCat), numpy.float32)

        # Test all the records to be sure that the measurement mechanism works for total flux
        # And that the area surrounding the footprint has the expected replacement pixels
        for i in range(len(measCat)):
            record = measCat[i]
            # First check to be sure that the flux measured by the plug-in is correct
            # This repeats the algorithm used by the NoiseReplacer, so not sure how useful it is
            foot = record.getFootprint()
            heavy = lsst.afw.detection.makeHeavyFootprint(foot, mi)
            noise = lsst.afw.detection.makeHeavyFootprint(foot, replaced)
            sumarray = numpy.ndarray((foot.getArea()), mi.getImage().getArray().dtype)
            lsst.afw.detection.flattenArray(foot, mi.getImage().getArray(), sumarray, mi.getImage().getXY0())
            sum = sumarray.sum(dtype=numpy.float64)
            count = foot.getArea()
            # get the values produced by the plugin 
            flux = record.get(fluxkey)
            fluxcount = record.get(fluxcountkey)
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
            arraysub = replaced.getImage().getArray()[ymin-y0:ymax-y0, xmin-x0:xmax-x0]
            bigflux = arraysub.sum(dtype=numpy.float64)
            
            # get the sum of the footprint area
            sumarray = numpy.ndarray((foot.getArea()), replaced.getImage().getArray().dtype)
            lsst.afw.detection.flattenArray(foot, replaced.getImage().getArray(), sumarray, replaced.getXY0())
            repflux = sumarray.sum(dtype=numpy.float64)

            newbackcount = (ymax-ymin)*(xmax-xmin) - count
            newback = bigflux - repflux
            back = record.get(backkey)
            backcount = record.get(backcountkey)


            # Test 1:  the area surrounding the object should be the same as what was measured
            #          during the actual run    
            self.assertEqual(backcount, newbackcount)
            self.assertEqual(back,newback)

            # Test 2:  the fill in of the area should be consistent with the noiseReplacer
            # mean=0.702828, std=27.4483
            popstd = math.sqrt((27.4483*27.4483)/len(sumarray))  
            normtest[i] = (sumarray.mean()-.702828)/popstd
        print normtest.max(), normtest.min(), normtest.mean(), normtest.std()

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
