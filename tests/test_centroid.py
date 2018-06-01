#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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
import unittest

import lsst.geom
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.utils.tests as utilsTests
import lsst.afw.detection.detectionLib as afwDetection
import lsst.afw.detection

try:
    type(verbose)
except NameError:
    display = False
    verbose = 0


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
        schema.getAliasMap().set("slot_Centroid", "test")
        centroider = alg(control, "test", schema)
        table = afwTable.SourceCatalog(schema)

        x0, y0 = 12345, 54321
        for imageFactory in (afwImage.MaskedImageF,):

            im = imageFactory(lsst.geom.ExtentI(100, 100))
            im.setXY0(lsst.geom.Point2I(x0, y0))

            # This fixed DoubleGaussianPsf replaces a computer generated one.
            # The values are not anything in particular, just a reasonable size.
            psf = lsst.afw.detection.GaussianPsf(15, 15, 3.0)
            exp = afwImage.makeExposure(im)
            exp.setPsf(psf)

            im.set(bkgd)
            x, y = 30, 20
            im.set(x, y, (1010,))

            source = table.addNew()
            spans = afwGeom.SpanSet(exp.getBBox(afwImage.LOCAL))
            foot = afwDetection.Footprint(spans)
            foot.addPeak(x + x0, y + y0, 1010)
            source.setFootprint(foot)
            centroider.measure(source, exp)

            self.assertFloatsAlmostEqual(x + x0, source.getX(), rtol=.00001)
            self.assertFloatsAlmostEqual(y + y0, source.getY(), rtol=.00001)
            self.assertFalse(source.get("test_flag"))

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
            self.assertFloatsAlmostEqual(x, source.getX(), rtol=.00001)
            self.assertFloatsAlmostEqual(y, source.getY(), rtol=.00001)
            self.assertFalse(source.get("test_flag"))

    def testNaiveMeasureCentroid(self):
        """Test that we can instantiate and play with NAIVE centroids"""
        bkgd = 10.0
        afwTable.SourceTable.makeMinimalSchema()
        control = measBase.NaiveCentroidControl()
        control.background = bkgd
        control.doFootprintCheck = False
        self.do_testAstrometry(measBase.NaiveCentroidAlgorithm, bkgd, control)

    def testSdssMeasureCentroid(self):
        """Test that we can instantiate and play with SDSS centroids"""
        control = measBase.SdssCentroidControl()
        control.doFootprintCheck = False
        self.do_testAstrometry(measBase.SdssCentroidAlgorithm, 10.0, control)


class SingleFrameMeasurementTaskTestCase(utilsTests.TestCase):
    """A test case for the SingleFrameMeasurementTask"""

    def mySetup(self):
        msConfig = measBase.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = ["base_SdssCentroid"]
        msConfig.slots.centroid = "base_SdssCentroid"
        msConfig.slots.shape = None
        msConfig.slots.psfShape = None
        msConfig.slots.apFlux = None
        msConfig.slots.modelFlux = None
        msConfig.slots.psfFlux = None
        msConfig.slots.instFlux = None
        msConfig.slots.calibFlux = None
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = measBase.SingleFrameMeasurementTask(schema, config=msConfig)
        measCat = afwTable.SourceCatalog(schema)

        source = measCat.addNew()
        spans = afwGeom.SpanSet(self.exp.getBBox(afwImage.LOCAL))
        fp = afwDetection.Footprint(spans)
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
        self.assertFloatsAlmostEqual(s.getCentroid(), lsst.geom.PointD(self.xcen, self.ycen), rtol=.01)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
