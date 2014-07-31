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
from lsst.meas.base.forcedCcd import ProcessForcedCcdTask, ProcessForcedCcdConfig
from lsst.daf.persistence.butler import *
import unittest
import lsst.utils.tests
import numpy

import lsst.afw.detection
import numpy

from lsst.meas.base.tests import *

numpy.random.seed(123)

#  Test SFM plugin, which is used to test that the plugin mechanism works correctly,
#  and that the noise replacement mechanism is working.  This plugin measures the total
#  flux for each source within its footprint, as well as the total flux in a box which
#  completely surrounds the object (providing an indication of any bad replacement nearby

class TestFlux(ForcedPlugin):

    executionOrder = 3.0

    def __init__(self, config, name, schemaMapper, flags=None, others=None, metadata=None):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.editOutputSchema()
        self.fluxKey = schema.addField("testFlux_flux", type=float, doc="sum of flux in object footprint",
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



ForcedPlugin.registry.register("testFlux", TestFlux)

DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")



def setUp():
    crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.005*lsst.afw.geom.degrees)
    catalog, bbox = MakeTestData.makeCatalog()
    exposure = MakeTestData.makeEmptyExposure(bbox, crval)
    MakeTestData.fillImages(catalog, exposure)
    catalog.writeFits(os.path.join(DATA_DIR, "truthcat.fits"))
    exposure.writeFits(os.path.join(DATA_DIR, "calexp.fits"))
    refCatalog = lsst.afw.table.SourceCatalog(catalog.getSchema())
    wcs = exposure.getWcs().clone()
    # shift the exposure to make a pretend reference exposure
    exposure.getWcs().shiftReferencePixel(lsst.afw.geom.Extent2D(1000,1000))
    refwcs = exposure.getWcs()
    exposure.writeFits(os.path.join(DATA_DIR, "refexp.fits"))
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
    catalog.writeFits(os.path.join(DATA_DIR, "refcat.fits"))
    refCat = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, 'refcat.fits'))
    centKey = lsst.afw.table.Point2DKey(catalog.getSchema().find("truth_x").key,
                                        catalog.getSchema().find("truth_y").key)
    coordKey = catalog.getSchema().find("coord").key 
    for rec in refCat:
        foot = rec.getFootprint()
        coord = rec.get(coordKey)
        cent = rec.get(centKey)
    crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 44.995*lsst.afw.geom.degrees)
    catalog, bbox = MakeTestData.makeCatalog()
    exposure = MakeTestData.makeEmptyExposure(bbox, crval)
    MakeTestData.fillImages(catalog, exposure)


def tearDown():
    os.unlink(os.path.join(DATA_DIR, "truthcat.fits"))
    os.unlink(os.path.join(DATA_DIR, "calexp.fits"))
    os.unlink(os.path.join(DATA_DIR, "refexp.fits"))
    os.unlink(os.path.join(DATA_DIR, "refcat.fits"))

#  Run the measurement (forcedCcd) task with the centroid peak plugin.
#  Measurement is done on a Ccd, which means that the coordinates will have
#  been transformed in generateSources.  Allow some tolerance for this
def testOnCcd():
    setUp()

    # get the reference catalog and reference exposure(might be a coadd)
    refCat = SourceCatalog.readFits(os.path.join(DATA_DIR, "refcat.fits"))
    path = os.path.join(DATA_DIR, 'refexp.fits')
    refexp = lsst.afw.image.ExposureF(path)
    refWcs = refexp.getWcs()

    # get the exposure to measure
    path = os.path.join(DATA_DIR, 'calexp.fits')
    exposure = lsst.afw.image.ExposureF(path)

    # configuration object for the forced task
    config = ForcedMeasurementConfig()
    config.plugins.names.clear()
    for plugin in ["centroid.peak", "testFlux"]:
        config.plugins.names.add(plugin)
    config.slots.centroid = "centroid.peak"
    config.slots.psfFlux = None
    config.slots.instFlux = "testFlux"
    config.slots.modelFlux = None
    config.slots.apFlux = None
    config.slots.shape = None

    # create and run the task
    task = ForcedMeasurementTask(refCat.getSchema(), config=config)
    result = task.run(exposure, refCat, refWcs)

    # Compare the results against the reference catalog
    sources = result.sources
    mismatches = 0
    testidkey = sources.getSchema().find("objectId").key
    truthFluxKey = refCat.getSchema().find("truth_flux").key
    for i, source in enumerate(sources):
        testFlux = sources[i].getInstFlux()
        truthFlux = refCat[i].get(truthFluxKey)
        parent = refCat[i].getParent()
        if parent==0:
            print source.getId(), truthFlux, testFlux
    tearDown()

if __name__ == "__main__":

    testOnCcd()
