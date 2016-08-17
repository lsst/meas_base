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

import unittest

import lsst.utils.tests
import lsst.meas.base.tests


class PixelFlagsTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.afw.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testNoFlags(self):
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(exposure, catalog)
        record = catalog[0]
        self.assertFalse(record.get("base_PixelFlags_flag"))
        self.assertFalse(record.get("base_PixelFlags_flag_edge"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolated"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturated"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_cr"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_bad"))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
