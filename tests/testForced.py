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
import unittest
import numpy

import lsst.utils.tests
import lsst.afw.detection
import lsst.afw.table

import lsst.meas.base.tests

numpy.random.seed(123)

# This task if just to define a ForcedMeasurementTask with no required plugins
class TestForcedMeasurementConfig(lsst.meas.base.ForcedMeasurementConfig):

    def setDefaults(self):
        super(lsst.meas.base.ForcedMeasurementConfig, self).setDefaults()
        self.plugins = ["base_PeakCentroid", "test_flux"]
        self.slots.centroid = None
        self.slots.shape = None
        self.slots.psfFlux = None
        self.slots.modelFlux = None
        self.slots.apFlux = None
        self.slots.instFlux = None

class TestForcedMeasurementTask(lsst.meas.base.ForcedMeasurementTask):

    ConfigClass = TestForcedMeasurementConfig


#  Test plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each source within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class TestFlux(lsst.meas.base.ForcedPlugin):

    def __init__(self, config, name, schemaMapper, flags=None, others=None, metadata=None):
        super(TestFlux, self).__init__(config, name, schemaMapper, flags,
                                       others, metadata)
        schema = schemaMapper.editOutputSchema()
        self.fluxKey = schema.addField("test_flux", type=float, doc="sum of flux in object footprint",
                                       units="dn")
        self.fluxCountKey = schema.addField("test_fluxcount", type=int,
                                            doc="number of pixels in object footprint", units="pixels^2")
        self.backKey = schema.addField("test_back", type=float, doc="avg of flux in background", units="dn")
        self.backCountKey = schema.addField("test_backcount", type=int,
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

lsst.meas.base.ForcedPlugin.registry.register("test_flux", TestFlux)

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")


class ForcedTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.005*lsst.afw.geom.degrees)
        catalog, bbox = lsst.meas.base.tests.MakeTestData.makeCatalog()
        exposure = lsst.meas.base.tests.MakeTestData.makeEmptyExposure(bbox, crval)
        lsst.meas.base.tests.MakeTestData.fillImages(catalog, exposure)
        catalog.writeFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "calexp-0A.fits"))
        refCatalog = lsst.afw.table.SourceCatalog(catalog.getSchema())
        wcs = exposure.getWcs().clone()
        exposure.getWcs().shiftReferencePixel(lsst.afw.geom.Extent2D(1000,1000))
        refwcs = exposure.getWcs()
        exposure.writeFits(os.path.join(DATA_DIR, "ref-0A.fits"))
        centKey = lsst.afw.table.Point2DKey(catalog.getSchema().find("truth_x").key,
                                            catalog.getSchema().find("truth_y").key)
        coordKey = catalog.getSchema().find("coord").key
        for rec in catalog:
            expRegion = exposure.getBBox(lsst.afw.image.PARENT)
            expRegion.shift(lsst.afw.geom.Extent2I(1000,1000))
            ref_footprint = rec.getFootprint().transform(wcs, refwcs, expRegion, True)
            peaks = rec.getFootprint().getPeaks()
            for i, peak in enumerate(peaks):
                newPeak = refwcs.skyToPixel(wcs.pixelToSky(lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())))
                ref_footprint.getPeaks().push_back(lsst.afw.detection.Peak(newPeak.getX(), newPeak.getY()))
            rec.setFootprint(ref_footprint)
            rfoot = rec.getFootprint()
            centroid = rec.get(centKey)
            sky = wcs.pixelToSky(centroid)
            rec.set(centKey, refwcs.skyToPixel(sky)) 
            rec.set(coordKey, sky) 
        catalog.writeFits(os.path.join(DATA_DIR, "refcat-0A.fits"))
        refCat = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, 'refcat-0A.fits'))
        centKey = lsst.afw.table.Point2DKey(catalog.getSchema().find("truth_x").key,
                                            catalog.getSchema().find("truth_y").key)
        coordKey = catalog.getSchema().find("coord").key 
        for rec in refCat:
            foot = rec.getFootprint()
            coord = rec.get(coordKey)
            cent = rec.get(centKey)
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 44.995*lsst.afw.geom.degrees)
        catalog, bbox = lsst.meas.base.tests.MakeTestData.makeCatalog()
        exposure = lsst.meas.base.tests.MakeTestData.makeEmptyExposure(bbox, crval)
        lsst.meas.base.tests.MakeTestData.fillImages(catalog, exposure)
        catalog.writeFits(os.path.join(DATA_DIR, "truthcat-0B.fits"))
        exposure.writeFits(os.path.join(DATA_DIR, "calexp-0B.fits"))
    

    def tearDown(self):
        os.unlink(os.path.join(DATA_DIR, "calexp-0B.fits"))
        os.unlink(os.path.join(DATA_DIR, "truthcat-0B.fits"))
        os.unlink(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "calexp-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "ref-0A.fits"))
        os.unlink(os.path.join(DATA_DIR, "refcat-0A.fits"))

    #  Run the measurement (forcedPhotCcd) task with the base_PeakCentroid plugin.
    #  Measurement is done on a Ccd, which means that the coordinates will have
    #  been transformed in generateSources.  Allow some tolerance for this
    def testOnCcd(self):

        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        exposure = lsst.afw.image.ExposureF(path)
        refCat = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "refcat-0A.fits"))
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
        task = TestForcedMeasurementTask(newRefCat.getSchema())
        result = task.run(exposure, list(newRefCat), refWcs)
        sources = result.sources
        mismatches = 0
	testidkey = sources.getSchema().find("objectId").key
	testFluxKey = sources.getSchema().find("test_flux").key
	truthFluxKey = refCat.getSchema().find("truth_flux").key
        for i, source in enumerate(sources):
            testFlux = sources[i].get(testFluxKey)
            truthFlux = refCat[i].get(truthFluxKey)
            parent = refCat[i].getParent()
            if parent==0 and abs((truthFlux-testFlux)/testFlux) > .03:
                mismatches += 1
        self.assertEqual(mismatches, 0) 

    #  Run the measurement (forcedPhotCcd) task with the base_PeakCentroid plugin.
    #  Measurement is done on the same exposure as was used for the reference catalog
    #  So in this case, the base_PeakCentroid should be exactly equal the peak[0] centroid
    #  Not that this does not have to be run on a Coadd for this test, just an exposure
    #  with the same Wcs as the refCat exposure
    def testOnSameWcs(self):

        refCat = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        refCat.table.defineCentroid("truth")
        refCat.table.defineShape("truth")
        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        exposure = lsst.afw.image.ExposureF(path)
        refWcs = exposure.getWcs()
        config = lsst.meas.base.ForcedMeasurementTask.ConfigClass()
        config.plugins = ["base_PeakCentroid"]
        config.slots.centroid = "base_PeakCentroid"
        config.slots.shape = None
        task = lsst.meas.base.ForcedMeasurementTask(config=config, refSchema=refCat.getSchema())
        result = task.run(exposure, refCat, refWcs)
        sources = result.sources
        mismatches = 0
        keyX = sources.getSchema().find("base_PeakCentroid_x").key
        keyY = sources.getSchema().find("base_PeakCentroid_y").key
        for source in sources:
            centroid = source.getFootprint().getPeaks()[0].getCentroid()
            peakcentroid = lsst.afw.geom.geomLib.Point2D(source.get(keyX), source.get(keyY))
            distance = numpy.sqrt(centroid.distanceSquared(peakcentroid))
            if distance > .00001:
                mismatches = mismatches + 1
                distX = centroid.getX()-peakcentroid.getX()
                distY = centroid.getY()-peakcentroid.getY()
                print "Mismatch: ", source.getId(), distance, distX, distY

        self.assertEqual(mismatches, 0)

    def testDefaultPlugins(self):

        refCat = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-0A.fits"))
        refCat.table.defineCentroid("truth")
        refCat.table.defineShape("truth")
        path = os.path.join(DATA_DIR, 'calexp-0A.fits')
        exposure = lsst.afw.image.ExposureF(path)
        refWcs = exposure.getWcs()
        config = lsst.meas.base.ForcedMeasurementTask.ConfigClass()
        task = lsst.meas.base.ForcedMeasurementTask(config=config, refSchema=refCat.getSchema())
        result = task.run(exposure, refCat, refWcs)
        sources = result.sources
        plugins = [
                 "base_GaussianFlux",
                 "base_NaiveFlux",
                 "base_PsfFlux",
                 "base_SincFlux",]
        for plugin in plugins:
            sources[0].get(plugin + "_flux")

        plugins = ["base_TransformedCentroid",]
        for plugin in plugins:
            sources[0].get(plugin + "_x")

        plugins = ["base_TransformedShape",]
        for plugin in plugins:
            sources[0].get(plugin + "_xx")



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
