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

import lsst.afw.geom
import lsst.meas.base.tests
import lsst.pex.exceptions
import lsst.utils.tests
import testLib  # noqa: F401 need this for SillyCentroid


class InputUtilitiesTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def testFlagAliases(self):
        """Test flag aliases are created correctly.

        In particular, we should get flag aliases to the slot centroid and
        shape algorithms when we initialize `GaussianFlux` (which uses both
        `SafeCentroidExtractor` and `SafeShapeExtractor`).
        """
        config = self.makeSingleFrameMeasurementConfig("base_GaussianFlux",
                                                       ["base_SdssCentroid", "base_SdssShape"])
        config.slots.centroid = "base_SdssCentroid"
        config.slots.shape = "base_SdssShape"
        task = self.makeSingleFrameMeasurementTask(config=config)
        # Test that the aliases resolve to the correct field.
        self.assertEqual(task.schema.find("base_GaussianFlux_flag_badCentroid").key,
                         task.schema.find("base_SdssCentroid_flag").key)
        self.assertEqual(task.schema.find("base_GaussianFlux_flag_badShape").key,
                         task.schema.find("base_SdssShape_flag").key)
        # Test that the aliases are direct links (i.e. they do not require
        # recursive expansion).
        self.assertEqual(task.schema.getAliasMap().get("base_GaussianFlux_flag_badCentroid"),
                         "base_SdssCentroid_flag")
        self.assertEqual(task.schema.getAliasMap().get("base_GaussianFlux_flag_badShape"),
                         "base_SdssShape_flag")

    def testCentroidFlagAliases(self):
        """Test aliases are correct when using multiple centroid algorithms.
        """
        config = self.makeSingleFrameMeasurementConfig("testLib_SillyCentroid", ["base_SdssCentroid"])
        config.slots.centroid = "base_SdssCentroid"
        config.slots.shape = None
        config.slots.psfShape = None
        task = self.makeSingleFrameMeasurementTask(config=config)
        # Test that the alias resolves to the correct field.
        self.assertEqual(task.schema.find("testLib_SillyCentroid_flag_badInitialCentroid").key,
                         task.schema.find("base_SdssCentroid_flag").key)
        # Test that the alias is a direct links (i.e. it do not require recursive expansion).
        self.assertEqual(task.schema.getAliasMap().get("testLib_SillyCentroid_flag_badInitialCentroid"),
                         "base_SdssCentroid_flag")
        # Test that there is no circular alias for the slot centroider itself.
        self.assertRaises(LookupError, task.schema.find, "base_SdssCentroid_flag_badInitialCentroid")

    def testUnmetCentroidDependency(self):
        """Test that an unmet centroid dependency raises.

        We should throw a `LogicError` when initializing an algorithm that
        requires a centroid without the centroid slot set.
        """
        config = self.makeSingleFrameMeasurementConfig("base_GaussianFlux",
                                                       ["base_SdssCentroid", "base_SdssShape"])
        config.slots.centroid = None
        config.slots.shape = "base_SdssShape"
        config.slots.psfShape = "base_SdssShape_psf"
        with self.assertRaises(lsst.pex.exceptions.LogicError):
            self.makeSingleFrameMeasurementTask(config=config)

    def testUnmetShapeDependency(self):
        """Test that an unmet shape dependency raises.

        Test that we throw a `LogicError` when initializing an algorithm that
        requires a shape without the shape slot set.
        """
        config = self.makeSingleFrameMeasurementConfig("base_GaussianFlux",
                                                       ["base_SdssCentroid", "base_SdssShape"])
        config.slots.centroid = "base_SdssCentroid"
        config.slots.shape = None
        config.slots.psfShape = None
        with self.assertRaises(lsst.pex.exceptions.LogicError):
            self.makeSingleFrameMeasurementTask(config=config)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
