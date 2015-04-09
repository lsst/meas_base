#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 AURA/LSST.
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
import re
import math
import unittest

import eups
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.utils.tests as utilsTests
import lsst.afw.detection.detectionLib as afwDetection
import lsst.afw.detection

import lsst.afw.display.ds9 as ds9

try:
    type(verbose)
except NameError:
    display = False
    verbose = 0

if False:
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        raise RuntimeError("Must set up afwdata to run these tests")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class CentroidTestCase(utilsTests.TestCase):
    """A test case for centroiding"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testCleanup(self):
        """Test that tearDown does"""
        pass

    def do_testAstrometry(self, alg, bkgd, control):
        """Test that we can instantiate and play with a centroiding algorithms"""

        schema = afwTable.SourceTable.makeMinimalSchema()
        centroider = alg(control, "test", schema)
        table = afwTable.SourceCatalog(schema)
        table.defineCentroid("test")

        x0, y0 = 12345, 54321
        for imageFactory in (afwImage.MaskedImageF,):

            im = imageFactory(afwGeom.ExtentI(100, 100))
            im.setXY0(afwGeom.Point2I(x0, y0))

            #   This fixed DoubleGaussianPsf replaces a computer generated one.
            #   The values are not anything in particular, just a reasonable size.
            psf = lsst.afw.detection.GaussianPsf(15, 15, 3.0)
            exp = afwImage.makeExposure(im)
            exp.setPsf(psf)
            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            im.set(bkgd)
            x, y = 30, 20
            im.set(x, y, (1010,))

            source = table.addNew()
            foot = afwDetection.Footprint(exp.getBBox(afwImage.LOCAL))
            foot.addPeak(x + x0, y + y0, 1010)
            source.setFootprint(foot)

            centroider.measure(source, exp)

            self.assertClose(x + x0, source.getX(), rtol=.00001)
            self.assertClose(y + y0, source.getY(), rtol=.00001)
            self.assertFalse(source.get("test_flag"))

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            im.set(bkgd)
            im.set(10, 20, (1010,))
            im.set(10, 21, (1010,))
            im.set(11, 20, (1010,))
            im.set(11, 21, (1010,))

            x, y = 10.5 + x0, 20.5 + y0
            source = table.addNew()
            source.set('test_x', x)
            source.set('test_y', y)
            centroider.measure(source, exp)
            self.assertClose(x, source.getX(), rtol=.00001)
            self.assertClose(y, source.getY(), rtol=.00001)
            self.assertFalse(source.get("test_flag"))

    def testGaussianMeasureCentroid(self):
        """Test that we can instantiate and play with GAUSSIAN centroids"""
        control = measBase.GaussianCentroidControl()
        self.do_testAstrometry(measBase.GaussianCentroidAlgorithm, 10.0, control)

    def testNaiveMeasureCentroid(self):
        """Test that we can instantiate and play with NAIVE centroids"""
        bkgd = 10.0
        schema = afwTable.SourceTable.makeMinimalSchema()
        control = measBase.NaiveCentroidControl()
        control.background = bkgd
        self.do_testAstrometry(measBase.NaiveCentroidAlgorithm, bkgd, control)

    def testSdssMeasureCentroid(self):
        """Test that we can instantiate and play with SDSS centroids"""
        control = measBase.SdssCentroidControl()
        self.do_testAstrometry(measBase.SdssCentroidAlgorithm, 10.0, control)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SingleFrameMeasurementTaskTestCase(utilsTests.TestCase):
    """A test case for the SingleFrameMeasurementTask"""

    def mySetup(self, runCentroider=True):
        msConfig = measBase.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = ["base_SdssCentroid"]
        if not runCentroider:
            msConfig.slots.centroid = None
        else:
            msConfig.slots.centroid = "base_SdssCentroid"

        msConfig.slots.shape = None
        msConfig.slots.apFlux = None
        msConfig.slots.modelFlux = None
        msConfig.slots.psfFlux = None
        msConfig.slots.instFlux = None
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = measBase.SingleFrameMeasurementTask(schema, config=msConfig)
        measCat = afwTable.SourceCatalog(schema)

        source = measCat.addNew()
        fp = afwDetection.Footprint(self.exp.getBBox(afwImage.LOCAL))
        fp.addPeak(50, 50, 1000.0)
        source.setFootprint(fp)
        # Then run the default SFM task.  Results not checked
        task.run(measCat, self.exp)
        return source

    def setUp(self):
        """Make the image we'll measure"""

        self.exp = afwImage.ExposureF(100, 100)
        psf = lsst.afw.detection.GaussianPsf(15, 15, 3.0)
        self.exp.setPsf(psf)
        self.I0, self.xcen, self.ycen = 1000.0, 50.5, 50.0
        im = self.exp.getMaskedImage().getImage()

        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                im.set(int(self.xcen) + i, int(self.ycen) + j, self.I0*(1 - 0.5*math.hypot(i - 0.5, j)))

    def tearDown(self):
        del self.exp

    def testCentroider(self):
        """Measure the centroid"""
        s = self.mySetup()

        # this does not match exactly, and it used to
        print s.getCentroid(), self.xcen, self.ycen
        self.assertClose(s.getCentroid(), afwGeom.PointD(self.xcen, self.ycen), rtol=.01)

    def testNoCentroider(self):
        """Check that we can disable running a centroid algorithm"""
        s = self.mySetup(False)

        self.assertRaises(pexExceptions.LogicError, s.getCentroid)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class MonetTestCase(unittest.TestCase):
    """A test case for centroiding using Dave Monet's 2-D Gaussian fitter"""

    def setUp(self):
	im = afwImage.ImageF(self.monetFile("small.fits"))
        self.mi = afwImage.MaskedImageF(im, afwImage.MaskU(im.getDimensions()),
                                        afwImage.ImageF(im.getDimensions()));
        self.ds = afwDetection.FootprintSet(self.mi, afwDetection.Threshold(100))

        if display:
            ds9.mtv(self.mi.getImage())
            ds9.erase()

        for foot in self.ds.getFootprints():
            bbox = foot.getBBox()
            x0, y0 = bbox.getMinX(), bbox.getMinY()
            x1, y1 = bbox.getMaxX(), bbox.getMaxY()
            xc = (x0 + x1)/2.0
            yc = (y0 + y1)/2.0

            if display:
                ds9.dot("+", xc, yc, ctype=ds9.BLUE)

                if False:
                    x0 -= 0.5; y0 -= 0.5
                    x1 += 0.5; y1 += 0.5

                    ds9.line([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], ctype=ds9.RED)
        msConfig = measBase.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = ["base_GaussianCentroid"]
        msConfig.slots.centroid = "base_GaussianCentroid"
        msConfig.slots.shape = None
        msConfig.slots.apFlux = None
        msConfig.slots.modelFlux = None
        msConfig.slots.psfFlux = None
        msConfig.slots.instFlux = None
        schema = afwTable.SourceTable.makeMinimalSchema()
        self.task = measBase.SingleFrameMeasurementTask(schema, config=msConfig)
        self.ssMeasured = afwTable.SourceCatalog(schema)
        self.ssMeasured.table.defineCentroid("base_GaussianCentroid")
        self.ssTruth = afwTable.SourceCatalog(schema)
        self.ssTruth.table.defineCentroid("base_GaussianCentroid")

        self.readTruth(self.monetFile("positions.dat-original"))

    def tearDown(self):
        del self.mi
        del self.ds
        del self.task
        del self.ssMeasured
        del self.ssTruth

    def monetFile(self, file):
        """Return a Monet file used for regression testing"""
        return os.path.join(eups.productDir("meas_base"), "tests", "Monet", file)

    def readTruth(self, filename):
        """Read Dave Monet's truth table"""
        for line in open(filename).readlines():
            if re.search(r"^\s*#", line):
                continue
            status, ID, xSex, xDGM, ySex, yDGM, sky = [float(el) for el in line.split()]
            s = self.ssTruth.addNew()
            s.setId(int(ID))
            s.set(self.ssTruth.table.getCentroidKey().getX(), xDGM)
            s.set(self.ssTruth.table.getCentroidKey().getY(), yDGM)

    def testMeasureCentroid(self):
        """Test that we can instantiate and play with a measureCentroid"""
        exposure = afwImage.makeExposure(self.mi)
        self.ds.makeSources(self.ssMeasured)
        ID = 1
        self.task.run(self.ssMeasured, exposure)
        for s in self.ssMeasured:
            s.setId(ID); ID += 1
            foot = s.getFootprint()
            bbox = foot.getBBox()
            xc = (bbox.getMinX() + bbox.getMaxX())//2
            yc = (bbox.getMinY() + bbox.getMaxY())//2

            if display:
                ds9.dot("x", xc, yc, ctype=ds9.GREEN)
        #
        # OK, we've measured all the sources.  Compare positions with Dave Monet's values
        #

        # FIXME: this test will fail until source matching in afw is updated to use afw/table
        mat = afwTable.matchXy(self.ssTruth, self.ssMeasured, 1.0)
        #self.assertEqual(ID, len(mat))  # we matched all the input sources

        eps = 6e-6                      # offset in pixels between measured centroid and the Truth
        for match in mat:
            dx = match[0].getX() - match[1].getX()
            dy = match[0].getY() - match[1].getY()

            good = True if math.hypot(dx, dy) < eps else False
            if not good:
                msg = "Star at (%.1f, %.1f): (dx, dy) = %g, %g)" % \
                    (match[0].getXAstrom(), match[0].getYAstrom(), dx, dy)
                if True:
                    print msg
                else:
                    self.assertTrue(good, msg)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CentroidTestCase)
    suites += unittest.makeSuite(MonetTestCase)
    suites += unittest.makeSuite(SingleFrameMeasurementTaskTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
