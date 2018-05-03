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

from contextlib import contextmanager
import unittest

import lsst.daf.base
import lsst.meas.base
import lsst.utils.tests

from lsst.meas.base.tests import (AlgorithmTestCase, CentroidTransformTestCase,
                                  SingleFramePluginTransformSetupHelper, ForcedPluginTransformSetupHelper)


@contextmanager
def onlyLogFatal(log):
    """
    For the duration of the context, only log FATAL errors.

    This is convenient when testing algorithms under failure conditions: we
    want to be able to check that they have set appropriate flags without
    spewing alarming & confusing error messages to the console.
    """
    oldLevel = log.getLevel()
    log.setLevel(log.FATAL)
    try:
        yield
    finally:
        log.setLevel(oldLevel)


class SingleFramePeakCentroidTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.afw.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(1000000.0, self.center)
        self.task = self.makeSingleFrameMeasurementTask("base_PeakCentroid")
        self.exposure, self.catalog = self.dataset.realize(10.0, self.task.schema, randomSeed=0)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset
        del self.task
        del self.exposure
        del self.catalog

    def testSingleFramePlugin(self):
        """Check that we recover the correct location of the centroid."""
        self.task.run(self.catalog, self.exposure)
        x = self.catalog[0].get("base_PeakCentroid_x")
        y = self.catalog[0].get("base_PeakCentroid_y")
        self.assertFalse(self.catalog[0].get("base_PeakCentroid_flag"))
        self.assertFloatsAlmostEqual(x, self.center.getX(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(y, self.center.getY(), atol=None, rtol=.02)

    def testFlags(self):
        """
        When it is impossible to measure the centroid -- in this case, because
        we have removed the Peaks from the SourceRecord -- the centroider
        should set a failure flag.
        """
        self.catalog[0].getFootprint().getPeaks().clear()
        # The decorator suppresses alarming but expected errors on the console.
        with onlyLogFatal(self.task.log):
            self.task.run(self.catalog, self.exposure)
        self.assertTrue(self.catalog[0].get("base_PeakCentroid_flag"))


class SingleFramePeakCentroidTransformTestCase(CentroidTransformTestCase,
                                               SingleFramePluginTransformSetupHelper,
                                               lsst.utils.tests.TestCase):

    class SingleFramePeakCentroidPluginFactory:
        """
        Helper class to sub in an empty PropertyList as the final argument to
        lsst.meas.base.SingleFramePeakCentroidPlugin.
        """

        def __call__(self, control, name, inputSchema):
            return lsst.meas.base.SingleFramePeakCentroidPlugin(control, name, inputSchema,
                                                                lsst.daf.base.PropertyList())
    controlClass = lsst.meas.base.SingleFramePeakCentroidConfig
    algorithmClass = SingleFramePeakCentroidPluginFactory()
    transformClass = lsst.meas.base.SimpleCentroidTransform
    flagNames = ()
    singleFramePlugins = ("base_PeakCentroid",)


class ForcedPeakCentroidTransformTestCase(CentroidTransformTestCase,
                                          ForcedPluginTransformSetupHelper,
                                          lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.ForcedPeakCentroidConfig
    algorithmClass = lsst.meas.base.ForcedPeakCentroidPlugin
    transformClass = lsst.meas.base.SimpleCentroidTransform
    flagNames = ()
    forcedPlugins = ("base_PeakCentroid",)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
