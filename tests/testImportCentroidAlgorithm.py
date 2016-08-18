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

from __future__ import division, absolute_import, print_function
import unittest

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.base
import lsst.utils.tests
import testLib

try:
    type(verbose)
except NameError:
    verbose = 0

lsst.meas.base.wrapSimpleAlgorithm(
    testLib.SillyCentroidAlgorithm,
    # algorithm name is specified manually here because testLib isn't in a normal Python package;
    # normally this is unnecessary
    name="testLib_SillyCentroid",
    Control=testLib.SillyCentroidControl,
    executionOrder=lsst.meas.base.BasePlugin.CENTROID_ORDER,
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class CentroidTestCase(lsst.utils.tests.TestCase):
    """A test case for centroiding."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testApplyCentroid(self):
        """Test that we can instantiate and play with SillyMeasureCentroid"""

        for imageFactory in (
            afwImage.MaskedImageF,
        ):
            im = imageFactory(afwGeom.ExtentI(100, 100))
            exp = afwImage.makeExposure(im)
            for offset in (0, 1, 2):
                control = testLib.SillyCentroidControl()
                control.param = offset
                x, y = 10, 20
                schema = afwTable.SourceTable.makeMinimalSchema()
                schema.addField("centroid_x", type=float)
                schema.addField("centroid_y", type=float)
                schema.getAliasMap().set("slot_Centroid", "centroid")
                plugin = testLib.SillyCentroidAlgorithm(control, "test", schema)
                measCat = afwTable.SourceCatalog(schema)
                source = measCat.makeRecord()
                source.set("centroid_x", x)
                source.set("centroid_y", y)
                plugin.measure(source, exp)
                self.assertEqual(x, source.get("test_x") - offset)
                self.assertEqual(y, source.get("test_y") - offset)

    def testMeasureCentroid(self):
        """Test that we can use our silly centroid through the usual Tasks."""
        testLib.SillyCentroidControl()
        x, y = 10, 20

        im = afwImage.MaskedImageF(afwGeom.ExtentI(512, 512))
        im.set(0)
        arr = im.getImage().getArray()
        arr[y, x] = 1000
        exp = afwImage.makeExposure(im)

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("flags_negative", type="Flag",
                        doc="set if source was detected as significantly negative")
        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        sfm_config.plugins = ["testLib_SillyCentroid"]
        sfm_config.plugins["testLib_SillyCentroid"].param = 5
        sfm_config.slots.centroid = "testLib_SillyCentroid"
        sfm_config.slots.shape = None
        sfm_config.slots.psfFlux = None
        sfm_config.slots.instFlux = None
        sfm_config.slots.calibFlux = None
        sfm_config.slots.apFlux = None
        sfm_config.slots.modelFlux = None
        sfm_config.doReplaceWithNoise = False
        task = lsst.meas.base.SingleFrameMeasurementTask(schema, config=sfm_config)
        measCat = afwTable.SourceCatalog(schema)
        measCat.defineCentroid("testLib_SillyCentroid")
        source = measCat.addNew()
        source.set("testLib_SillyCentroid_x", x)
        source.set("testLib_SillyCentroid_y", y)
        source.set("parent", 0)
        source.set("flags_negative", False)

        # now run the SFM task with the test plugin
        task.run(measCat, exp)
        self.assertEqual(len(measCat), 1)
        self.assertEqual(measCat[0].getY(), y + 5)
        self.assertEqual(measCat[0].getX(), x + 5)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
