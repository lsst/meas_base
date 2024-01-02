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

import numpy as np

import lsst.geom
import lsst.utils.tests
import lsst.meas.base.tests


class PixelFlagsTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def testNoFlags(self):
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]
        self.assertFalse(record.get("base_PixelFlags_flag"))
        self.assertFalse(record.get("base_PixelFlags_flag_edge"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolated"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturated"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_saturatedCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_cr"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenterAll"))
        self.assertFalse(record.get("base_PixelFlags_flag_bad"))
        self.assertFalse(record.get("base_PixelFlags_flag_badCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_badCenterAll"))

    def testSomeFlags(self):
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        # one cr pixel outside the center
        cosmicray = exposure.mask.getPlaneBitMask("CR")
        x = round(self.center.x)
        y = round(self.center.y)
        exposure.mask[x+3, y+4] |= cosmicray
        # one interpolated pixel near the center
        interpolated = exposure.mask.getPlaneBitMask("INTRP")
        exposure.mask[self.center] |= interpolated
        # all pixels in the center are bad
        bad = exposure.mask.getPlaneBitMask("BAD")
        exposure.mask[x-1:x+2, y-1:y+2] |= bad
        task.run(catalog, exposure)
        record = catalog[0]

        self.assertTrue(record.get("base_PixelFlags_flag_cr"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_crCenterAll"))
        self.assertTrue(record.get("base_PixelFlags_flag_interpolated"))
        self.assertTrue(record.get("base_PixelFlags_flag_interpolatedCenter"))
        self.assertFalse(record.get("base_PixelFlags_flag_interpolatedCenterAll"))
        self.assertTrue(record.get("base_PixelFlags_flag_bad"))
        self.assertTrue(record.get("base_PixelFlags_flag_badCenter"))
        self.assertTrue(record.get("base_PixelFlags_flag_badCenterAll"))

    def testOffimageNonfinite(self):
        """Test that flag_offimage and flag_edge are set correctly for
        nonfinite sources.
        """
        # These three will be explicitly set to have a non-finite centroid; we
        # cannot add sources off the image in TestDataset.
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 20))
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 100))
        self.dataset.addSource(100000, lsst.geom.Point2D(100, 30))
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        catalog[1]["slot_Centroid_x"] = np.inf
        catalog[2]["slot_Centroid_x"] = -np.inf
        catalog[3]["slot_Centroid_y"] = np.nan
        task.run(catalog, exposure)
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_offimage"],
                                      np.array([False, True, True, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edge"],
                                      np.array([False, True, True, True]))

    def testOffimage(self):
        """Test that sources at the boundary of the image get flag_offimage
        and flag_edge set appropriately.
        """
        # These four will be explicitly set to values at the boundary; we
        # cannot add sources off the image in TestDataset.
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 20))
        self.dataset.addSource(100000, lsst.geom.Point2D(20, 100))
        self.dataset.addSource(100000, lsst.geom.Point2D(80, 100))
        self.dataset.addSource(100000, lsst.geom.Point2D(80, 30))
        task = self.makeSingleFrameMeasurementTask("base_PixelFlags")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        catalog[1]["slot_Centroid_x"] = -20.5  # on image
        catalog[2]["slot_Centroid_x"] = -20.51  # off image
        catalog[3]["slot_Centroid_y"] = 129.4  # on image
        catalog[4]["slot_Centroid_y"] = 129.50  # off image
        task.run(catalog, exposure)
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_offimage"],
                                      np.array([False, False, True, False, True]))
        np.testing.assert_array_equal(catalog["base_PixelFlags_flag_edge"],
                                      np.array([False, False, True, False, True]))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
