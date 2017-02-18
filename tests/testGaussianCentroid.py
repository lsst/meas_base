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

from __future__ import absolute_import, division, print_function
import unittest

import lsst.utils.tests
import lsst.meas.base.tests
from lsst.meas.base.tests import (AlgorithmTestCase, CentroidTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)


class GaussianCentroidTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

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

    def testSingleFramePlugin(self):
        task = self.makeSingleFrameMeasurementTask("base_GaussianCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(catalog, exposure)
        record = catalog[0]
        x = record.get("base_GaussianCentroid_x")
        y = record.get("base_GaussianCentroid_y")
        self.assertFalse(record.get("base_GaussianCentroid_flag"))
        self.assertFalse(record.get("base_GaussianCentroid_flag_noPeak"))
        self.assertClose(x, self.center.getX(), atol=None, rtol=.01)
        self.assertClose(y, self.center.getY(), atol=None, rtol=.01)


class GaussianCentroidTransformTestCase(CentroidTransformTestCase, SingleFramePluginTransformSetupHelper,
                                        lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.GaussianCentroidControl
    algorithmClass = lsst.meas.base.GaussianCentroidAlgorithm
    transformClass = lsst.meas.base.GaussianCentroidTransform
    flagNames = ('flag', 'flag_noPeak')
    singleFramePlugins = ('base_GaussianCentroid',)
    forcedPlugins = ('base_GaussianCentroid',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
