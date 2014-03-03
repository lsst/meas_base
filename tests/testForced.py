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
from lsst.afw.table import IdFactory,Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.base import *
from lsst.meas.base.forcedImage import *
from lsst.meas.base.forcedCcd import ForcedCcdMeasurementTask, ForcedCcdMeasurementConfig
from lsst.daf.persistence.butler import *
import unittest
import lsst.utils.tests
import numpy

import lsst.afw.detection
import numpy

from lsst.meas.base.tests import *

numpy.random.seed(123)

class TestForcedCcdMeasurementConfig(ForcedCcdMeasurementConfig):

    def setDefaults(self):
        ForcedCcdMeasurementConfig.setDefaults(self)
        self.plugins = ["centroid.peak", "test.flux"]
        self.slots.centroid = "centroid.peak"
        self.slots.shape = None
        self.slots.psfFlux = None 
        self.slots.modelFlux = None
        self.slots.apFlux = None
        self.slots.instFlux = None

class TestForcedCcdMeasurementTask(ForcedCcdMeasurementTask):

    ConfigClass = TestForcedCcdMeasurementConfig
    # Since this is a test and we are not building a real output catalog
    # we can use create Ids without using camera dependent code.
    def makeIdFactory(self, butler):
        return IdFactory.makeSimple()

class TestFluxConfig(ForcedPluginConfig):
    pass

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each source within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class TestFlux(ForcedPlugin):
    ConfigClass = TestFluxConfig

    def __init__(self, config, name, schemaMapper, flags=None, others=None, metadata=None):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.editOutputSchema()
        self.fluxKey = schema.addField("test.flux", type=float, doc="sum of flux in object footprint",
                                       units="dn")
        self.fluxCountKey = schema.addField("test.fluxcount", type=int,
                                            doc="number of pixels in object footprint", units="pixels^2")
        self.backKey = schema.addField("test.back", type=float, doc="avg of flux in background", units="dn")
        self.backCountKey = schema.addField("test.backcount", type=int,
                                            doc="number of pixels in surrounding background",
                                            units="pixels^2")

    def measure(self, measRecord, exposure, refRecord, referenceWcs):
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



ForcedPlugin.registry.register("test.flux", TestFlux)

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests", "data")


class ForcedTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.005*lsst.afw.geom.degrees)
        catalog, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox, crval)
        MakeTestData.fillImages(catalog, exposure)
        catalog.writeFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "calexp-0A.fits"))
        refcatalog = lsst.afw.table.SourceCatalog(catalog.getSchema())
        wcs = exposure.getWcs().clone()
        exposure.getWcs().shiftReferencePixel(lsst.afw.geom.Extent2D(1000,1000))
        refwcs = exposure.getWcs()
        exposure.writeFits(os.path.join(DATA_DIR, "ref-0A.fits"))
        centkey = catalog.getSchema().find("truth.centroid").key 
        coordkey = catalog.getSchema().find("coord").key 
        for rec in catalog:
            expregion = exposure.getBBox(lsst.afw.image.PARENT)
            expregion.shift(lsst.afw.geom.Extent2I(1000,1000))
            ref_footprint = rec.getFootprint().transform(wcs, refwcs, expregion, True)
            peaks = rec.getFootprint().getPeaks()
            for i, peak in enumerate(peaks):
                newpeak = refwcs.skyToPixel(wcs.pixelToSky(lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())))
                ref_footprint.getPeaks().push_back(lsst.afw.detection.Peak(newpeak.getX(), newpeak.getY()))
            rec.setFootprint(ref_footprint)
            rfoot = rec.getFootprint()
            centroid = rec.get(centkey)
            sky = wcs.pixelToSky(centroid)
            rec.set(centkey, refwcs.skyToPixel(sky)) 
            rec.set(coordkey, sky) 
        catalog.writeFits(os.path.join(DATA_DIR, "refcat-0A.fits"))
        refCat = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, 'refcat-0A.fits'))
        centkey = catalog.getSchema().find("truth.centroid").key 
        coordkey = catalog.getSchema().find("coord").key 
        for rec in refCat:
            foot = rec.getFootprint()
            coord = rec.get(coordkey)
            cent = rec.get(centkey)
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 44.995*lsst.afw.geom.degrees)
        catalog, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox, crval)
        MakeTestData.fillImages(catalog, exposure)
        catalog.writeFits(os.path.join(DATA_DIR, "truthcat-0B.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "calexp-0B.fits"))
    

    def tearDown(self):
        os.unlink(os.path.join(DATA_DIR, "calexp-0B.fits"))
        os.unlink(os.path.join(DATA_DIR, "truthcat-0B.fits"))
        os.unlink(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "calexp-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "ref-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "refcat-0A.fits"))

    #  Run the measurement (forcedCcd) task with the centroid peak plugin.
    #  Measurement is done on a Ccd, which means that the coordinates will have
    #  been transformed in generateSources.  Allow some tolerance for this
    def testOnCcd(self):

        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        exposure = lsst.afw.image.ExposureF(path)
        refCat = SourceCatalog.readFits(os.path.join(DATA_DIR, "refcat-0A.fits"))
        path = os.path.join(DATA_DIR, 'ref-0A.fits')
        ref = lsst.afw.image.ExposureF(path)
        refWcs = ref.getWcs()
        rec = refCat[1]
        foot = rec.getFootprint()
        mapper = lsst.afw.table.SchemaMapper(refCat.getSchema())
        minimalSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        mapper.addMinimalSchema(minimalSchema)
        newRefCat = lsst.afw.table.SourceCatalog(mapper.getOutputSchema()) 
        newRefCat.extend(refCat, mapper=mapper)
        # Read a catalog which should be the same as the catalog of processCcd
        # prior to measurement.  Create an empty catalog with the same schema
        # plus the schema items for the SFM task, then transfer the existing data
        # to the new catalog
        srccat = SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        schema = srccat.getSchema()
        flags = MeasurementDataFlags()
        
        task = TestForcedCcdMeasurementTask(refSchema=newRefCat.getSchema())

        # Then run the default SFM task.  Results not checked
        result = task.forcedMeasure(exposure, list(newRefCat), refWcs)
        sources = result.sources
        mismatches = 0
	testidkey = sources.getSchema().find("objectId").key
	truthidkey = refCat.getSchema().find("id").key
	testfluxkey = sources.getSchema().find("test.flux").key
	truthfluxkey = refCat.getSchema().find("truth.flux").key
        for i, source in enumerate(sources):
            testflux = sources[i].get(testfluxkey)
            truthflux = refCat[i].get(truthfluxkey)
            parent = refCat[i].getParent()
            if parent==0 and abs((truthflux-testflux)/testflux) > .03:
                mismatches += 1
        self.assertEqual(mismatches, 0) 

    #  Run the measurement (forcedCcd) task with the centroid peak plugin.
    #  Measurement is done on the same exposure as was used for the reference catalog
    #  So in this case, the centroid.peak should be exactly equal the peak[0] centroid
    #  Not that this does not have to be run on a Coadd for this test, just an exposure
    #  with the same Wcs as the refCat exposure
    def testOnSameWcs(self):

        refCat = SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        
        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        exposure = lsst.afw.image.ExposureF(path)
        refWcs = exposure.getWcs()

        task = TestForcedCcdMeasurementTask(butler=None, refSchema=refCat.getSchema())
        result = task.forcedMeasure(exposure, refCat, refWcs) #butler
        sources = result.sources
        mismatches = 0
        key = sources.getSchema().find("centroid.peak").key
        for source in sources:
            centroid = source.getFootprint().getPeaks()[0].getCentroid()
            peakcentroid = source.get(key)
            distance = math.sqrt(centroid.distanceSquared(peakcentroid))
            if distance > .00001:
                mismatches = mismatches + 1
                distX = centroid.getX()-peakcentroid.getX()
                distY = centroid.getY()-peakcentroid.getY()
                print "Mismatch: ", source.getId(), distance, distX, distY

        self.assertEqual(mismatches, 0) 

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ForcedTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
