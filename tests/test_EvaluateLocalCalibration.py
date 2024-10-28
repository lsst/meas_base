# This file is part of ap_association.
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

import logging
import numpy as np
import unittest

import lsst.geom
import lsst.utils.tests
import lsst.meas.base.tests


class TestLocalPhotoCalibration(lsst.meas.base.tests.AlgorithmTestCase,
                                lsst.utils.tests.TestCase):

    def setUp(self):
        self.pluginName = "base_LocalPhotoCalib"
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testPhotoCalib(self):
        task = self.makeSingleFrameMeasurementTask(self.pluginName)
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]

        calib = exposure.getPhotoCalib().getLocalCalibration(self.center)
        calibErr = exposure.getPhotoCalib().getCalibrationErr()
        self.assertEqual(record.get(self.pluginName), calib)
        self.assertEqual(record.get(self.pluginName + "Err"), calibErr)

    def testPhotoCalibIsNone(self):
        with self.assertNoLogs(level=logging.WARNING) and self.assertLogs(level=logging.DEBUG):
            task = self.makeSingleFrameMeasurementTask(self.pluginName)
            exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
            exposure.setPhotoCalib(None)
            task.run(catalog, exposure)
            record = catalog[0]

            self.assertTrue(np.isnan(record.get(self.pluginName)))
            self.assertTrue(np.isnan(record.get(self.pluginName + "Err")))
            self.assertTrue(record.get(self.pluginName + "_flag"))


class TestLocalWcs(lsst.meas.base.tests.AlgorithmTestCase,
                   lsst.utils.tests.TestCase):

    def setUp(self):
        self.pluginName = "base_LocalWcs"
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testMatrix(self):
        task = self.makeSingleFrameMeasurementTask(self.pluginName)
        exposure, catalog = self.dataset.realize(10.0,
                                                 task.schema,
                                                 randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]

        # Get CdMatrix from afw. Convert from degrees to radians.
        trueCdMatrix = np.radians(exposure.getWcs().getCdMatrix())
        self.assertAlmostEqual(record.get(self.pluginName + "_CDMatrix_1_1"),
                               trueCdMatrix[0, 0])
        self.assertAlmostEqual(record.get(self.pluginName + "_CDMatrix_2_1"),
                               trueCdMatrix[1, 0])
        self.assertAlmostEqual(record.get(self.pluginName + "_CDMatrix_1_2"),
                               trueCdMatrix[0, 1])
        self.assertAlmostEqual(record.get(self.pluginName + "_CDMatrix_2_2"),
                               trueCdMatrix[1, 1])
        self.assertAlmostEqual(
            exposure.getWcs().getPixelScale(record.getCentroid()).asRadians(),
            np.sqrt(np.fabs(record.get(self.pluginName + "_CDMatrix_1_1")
                            * record.get(self.pluginName + "_CDMatrix_2_2")
                            - record.get(self.pluginName + "_CDMatrix_2_1")
                            * record.get(self.pluginName + "_CDMatrix_1_2"))))

    def testWcsIsNone(self):
        with self.assertNoLogs(level=logging.WARNING) and self.assertLogs(level=logging.DEBUG):
            task = self.makeSingleFrameMeasurementTask(self.pluginName)
            exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
            exposure.setWcs(None)
            task.run(catalog, exposure)
            record = catalog[0]

            self.assertTrue(np.isnan(record.get(self.pluginName + "_CDMatrix_1_1")))
            self.assertTrue(np.isnan(record.get(self.pluginName + "_CDMatrix_2_1")))
            self.assertTrue(np.isnan(record.get(self.pluginName + "_CDMatrix_1_2")))
            self.assertTrue(np.isnan(record.get(self.pluginName + "_CDMatrix_2_2")))
            self.assertTrue(record.get(self.pluginName + "_flag"))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
