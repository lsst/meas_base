# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import lsst.geom
from lsst.meas.base.tests import (AlgorithmTestCase, CentroidTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)
import lsst.utils.tests


class NaiveCentroidTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testSingleFramePlugin(self):
        task = self.makeSingleFrameMeasurementTask("base_NaiveCentroid")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]
        x = record.get("base_NaiveCentroid_x")
        y = record.get("base_NaiveCentroid_y")
        self.assertFalse(record.get("base_NaiveCentroid_flag"))
        self.assertFloatsAlmostEqual(x, self.center.getX(), atol=2.)
        self.assertFloatsAlmostEqual(y, self.center.getY(), atol=2.)


class NaiveCentroidTransformTestCase(CentroidTransformTestCase,
                                     SingleFramePluginTransformSetupHelper,
                                     lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.NaiveCentroidControl
    algorithmClass = lsst.meas.base.NaiveCentroidAlgorithm
    transformClass = lsst.meas.base.NaiveCentroidTransform
    flagNames = ('flag', 'flag_noCounts', 'flag_edge')
    singleFramePlugins = ('base_NaiveCentroid',)
    forcedPlugins = ('base_NaiveCentroid',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
