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

import numpy as np
import pandas as pd
import unittest

from lsst.meas.base import (
    DiaObjectCalculationTask,
    DiaObjectCalculationConfig,
    DiaObjectCalculationPlugin)
from lsst.meas.base.pluginRegistry import register
import lsst.utils.tests


@register("testCount")
class CountDiaPlugin(DiaObjectCalculationPlugin):
    """Simple mean function.
    """
    outputCols = ["count"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaObjectId,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """
        """
        diaObjects.at[diaObjectId, "count"] = len(diaSources["psfFlux"])


@register("testDiaPlugin")
class DiaPlugin(DiaObjectCalculationPlugin):
    """Simple mean function.
    """
    outputCols = ["MeanFlux", "StdFlux"]

    plugType = "multi"

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_CATALOGCALCULATION

    def calculate(self,
                  diaObjects,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        """
        """
        diaObjects.loc[:, "%sMeanFlux" % band] = \
            filterDiaSources.psfFlux.agg("mean")
        diaObjects.loc[:, "%sStdFlux" % band] = \
            filterDiaSources.psfFlux.agg("std")


@register("testDependentDiaPlugin")
class DependentDiaPlugin(DiaObjectCalculationPlugin):
    """Simple calculation using the previously calculated mean.
    """
    inputCols = ["MeanFlux"]
    outputCols = ["ChiFlux"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self,
                  diaObjects,
                  diaObjectId,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        diaObjects.at[diaObjectId, "%sChiFlux" % band] = np.sum(
            ((filterDiaSources["psfFlux"]
              - diaObjects.at[diaObjectId, "%sMeanFlux" % band])
             / filterDiaSources["psfFluxErr"]) ** 2)


@register("testCollidingDiaPlugin")
class CollidingDiaPlugin(DiaObjectCalculationPlugin):
    """Simple calculation using the previously calculated mean.
    """
    outputCols = ["MeanFlux"]

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_MOMENTS_CALCULATED

    def calculate(self,
                  diaObjects,
                  diaObjectId,
                  diaSources,
                  filterDiaSources,
                  band,
                  **kwargs):
        diaObjects.at[diaObjectId, "%sMeanFlux" % band] = 0.0


class TestDiaCalcluation(unittest.TestCase):

    def setUp(self):
        # Create diaObjects
        self.newDiaObjectId = 13
        self.diaObjects = pd.DataFrame(
            data=[{"diaObjectId": objId}
                  for objId in [0, 1, 2, 3, 4, 5, self.newDiaObjectId]])

        # Create diaSources from "previous runs" and newly created ones.
        diaSources = [{"diaSourceId": objId, "diaObjectId": objId,
                       "psfFlux": 0., "psfFluxErr": 1.,
                       "scienceFlux": 0., "scienceFluxErr": 1.,
                       "midpointMjdTai": 0, "band": "g"}
                      for objId in range(5)]
        diaSources.extend([{"diaSourceId": 5 + objId, "diaObjectId": objId,
                            "psfFlux": 0., "psfFluxErr": 1.,
                            "scienceFlux": 0., "scienceFluxErr": 1.,
                            "midpointMjdTai": 0, "band": "r"}
                           for objId in range(5)])
        diaSources.extend([{"diaSourceId": 10, "diaObjectId": 0,
                            "psfFlux": 1., "psfFluxErr": 1.,
                            "scienceFlux": 0., "scienceFluxErr": 0.,
                            "midpointMjdTai": 0, "band": "g"},
                           {"diaSourceId": 11, "diaObjectId": 1,
                            "psfFlux": 1., "psfFluxErr": 1.,
                            "scienceFlux": 0., "scienceFluxErr": 0.,
                            "midpointMjdTai": 0, "band": "g"},
                           {"diaSourceId": 12, "diaObjectId": 2,
                            "psfFlux": np.nan, "psfFluxErr": 1.,
                            "scienceFlux": 0., "scienceFluxErr": 0.,
                            "midpointMjdTai": 0, "band": "g"},
                           {"diaSourceId": self.newDiaObjectId,
                            "diaObjectId": self.newDiaObjectId,
                            "psfFlux": 1., "psfFluxErr": 1.,
                            "scienceFlux": 0., "scienceFluxErr": 0.,
                            "midpointMjdTai": 0, "band": "g"}])
        self.diaSources = pd.DataFrame(data=diaSources)

        self.updatedDiaObjectIds = np.array([0, 1, 2, self.newDiaObjectId],
                                            dtype=np.int64)

        conf = DiaObjectCalculationConfig()
        conf.plugins = ["testDiaPlugin",
                        "testDependentDiaPlugin"]
        self.diaObjCalTask = DiaObjectCalculationTask(config=conf)

    def testRun(self):
        """Test the run method and that diaObjects are updated correctly.
        """
        results = self.diaObjCalTask.run(self.diaObjects,
                                         self.diaSources,
                                         self.updatedDiaObjectIds,
                                         ["g"])
        diaObjectCat = results.diaObjectCat
        updatedDiaObjects = results.updatedDiaObjects
        updatedDiaObjects.set_index("diaObjectId", inplace=True)
        # Test the lengths of the output dataframes.
        self.assertEqual(len(diaObjectCat), len(self.diaObjects))
        self.assertEqual(len(updatedDiaObjects),
                         len(self.updatedDiaObjectIds))

        # Test values stored computed in the task.
        for objId, diaObject in updatedDiaObjects.iterrows():
            if objId == self.newDiaObjectId:
                self.assertEqual(diaObject["gMeanFlux"], 1.)
                self.assertTrue(np.isnan(diaObject["gStdFlux"]))
                self.assertAlmostEqual(diaObject["gChiFlux"], 0.0)
            elif objId == 2:
                self.assertAlmostEqual(diaObject["gMeanFlux"], 0.0)
                self.assertTrue(np.isnan(diaObject["gStdFlux"]))
                self.assertAlmostEqual(diaObject["gChiFlux"], 0.0)
            else:
                self.assertAlmostEqual(diaObject["gMeanFlux"], 0.5)
                self.assertAlmostEqual(diaObject["gStdFlux"],
                                       0.7071067811865476)
                self.assertAlmostEqual(diaObject["gChiFlux"], 0.5)

    def testRunUnindexed(self):
        """Test inputing un-indexed catalogs.
        """
        unindexedDiaSources = pd.DataFrame(data=[
            {"diaSourceId": objId, "diaObjectId": 0,
             "psfFlux": 0., "psfFluxErr": 1.,
             "scienceFlux": 0., "scienceFluxErr": 1.,
             "midpointMjdTai": 0, "band": "g"}
            for objId in range(1000)])
        unindexedDiaSources = pd.concat(
            (
                unindexedDiaSources,
                pd.DataFrame(
                    data=[
                        {
                            "diaSourceId": objId + 1000,
                            "diaObjectId": 0,
                            "psfFlux": 0., "psfFluxErr": 1.,
                            "scienceFlux": 0., "scienceFluxErr": 1.,
                            "midpointMjdTai": 0, "band": "g",
                        }
                        for objId in range(10)
                    ]
                )
            )
        )

        conf = DiaObjectCalculationConfig()
        conf.plugins = ["testCount"]
        diaObjectCalTask = DiaObjectCalculationTask(config=conf)
        self.diaObjects.reset_index()
        results = diaObjectCalTask.run(self.diaObjects,
                                       unindexedDiaSources,
                                       np.array([0], dtype=np.int64),
                                       ["g"])
        updatedDiaObjects = results.updatedDiaObjects
        self.assertEqual(updatedDiaObjects.at[0, "count"],
                         len(unindexedDiaSources))

    def testConflictingPlugins(self):
        """Test that code properly exits upon plugin collision.
        """
        with self.assertRaises(ValueError):
            conf = DiaObjectCalculationConfig()
            conf.plugins = ["testDependentDiaPlugin"]
            DiaObjectCalculationTask(config=conf)

        with self.assertRaises(ValueError):
            conf = DiaObjectCalculationConfig()
            conf.plugins = ["testDiaPlugin",
                            "testCollidingDiaPlugin",
                            "testDependentDiaPlugin"]
            DiaObjectCalculationTask(config=conf)

        # Test that ordering in the config does not matter and dependent
        # plugin is instantiated after independent plugin. Would raise
        # ValueError on failure.
        conf = DiaObjectCalculationConfig()
        conf.plugins = ["testDependentDiaPlugin",
                        "testDiaPlugin"]
        DiaObjectCalculationTask(config=conf)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
