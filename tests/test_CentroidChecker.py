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

import lsst.utils.tests
import lsst.geom
import lsst.meas.base
import lsst.meas.base.tests
from lsst.meas.base.tests import (AlgorithmTestCase)
import lsst.pex.config
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base import FlagDefinitionList, FlagHandler


class CentroiderConfig(SingleFramePluginConfig):

    moveX = lsst.pex.config.Field(dtype=int, default=0,
                                  doc="amount to re-position in X")
    moveY = lsst.pex.config.Field(dtype=int, default=0,
                                  doc="amount to re-position in Y")
    dist = lsst.pex.config.Field(dtype=int, default=0,
                                 doc="distance to allow centroid to be off")
    setErrors = lsst.pex.config.Field(dtype=bool, default=False,
                                      doc="set errors on measurement to errX, errY")
    errX = lsst.pex.config.Field(dtype=float, default=0,
                                 doc="uncertainty on X measurement")
    errY = lsst.pex.config.Field(dtype=float, default=0,
                                 doc="uncertainty on X measurement")


@register("test_Centroider")
class Centroider(SingleFramePlugin):
    """Sample Python measurement plugin.

    The flag handler for this plugin is created during construction, and is
    called using the method `fail`.  All plugins are required to implement
    this method, which is used to set the flags in the output source record if
    an error occurs.
    """
    ConfigClass = CentroiderConfig
    # Class variables ErrEnum and FLAGDEFS are added by the decorator

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)

        flagDefs = FlagDefinitionList()
        flagDefs.add("flag", "General Failure error")
        flagDefs.add("test_flag", "second flag")
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)

        if self.config.setErrors:
            uncertainty = lsst.meas.base.UncertaintyEnum.SIGMA_ONLY
        else:
            uncertainty = lsst.meas.base.UncertaintyEnum.NO_UNCERTAINTY

        self.centroidKey = lsst.meas.base.CentroidResultKey.addFields(schema, name, name, uncertainty)

        if self.config.dist is None:
            self.centroidChecker = lsst.meas.base.CentroidChecker(schema, name)
        else:
            self.centroidChecker = lsst.meas.base.CentroidChecker(schema, name, True, self.config.dist)

    def measure(self, measRecord, exposure):
        """This measure routine moves the centroid by design to create an error.
        """
        measRecord.set(self.centroidKey.getX(), measRecord.getX() + self.config.moveX)
        measRecord.set(self.centroidKey.getY(), measRecord.getY() + self.config.moveY)
        if self.centroidKey.getCentroidErr().isValid():
            err = measRecord.get(self.centroidKey.getCentroidErr())
            err[0][0] = self.config.errX
            err[1][1] = self.config.errY
            measRecord.set(self.centroidKey.getCentroidErr(), err)
        self.centroidChecker(measRecord)

    def fail(self, measRecord, error=None):
        """Respond to measurement failures.

        This routine responds to the standard failure call in baseMeasurement
        If the exception is a MeasurementError, the error will be passed to
        the fail method by the MeasurementFramework.
        """
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)


class CentroidCheckerTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    # Setup a configuration and datasource to be used by the plugin tests
    def setUp(self):
        self.algName = "test_Centroider"
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(bbox)
        self.dataset.addSource(instFlux=1E5, centroid=lsst.geom.Point2D(25, 26))

    def tearDown(self):
        del self.dataset

    def testNoError(self):
        """Test that the ``resetToPeak`` flag is not set when no error seen.
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeSingleFrameMeasurementConfig(plugin=self.algName)
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=0)
        task.run(cat, exposure)
        source = cat[0]
        self.assertFalse(source.get("test_Centroider_flag"))
        self.assertFalse(source.get("test_Centroider_flag_resetToPeak"))

    def testCheckErrors(self):
        """Test that centroids with invalid (NaN) errors are flagged.
        """
        def runMeasurement(errX, errY):
            schema = self.dataset.makeMinimalSchema()
            config = self.makeSingleFrameMeasurementConfig(plugin=self.algName)
            config.plugins[self.algName].setErrors = True
            config.plugins[self.algName].errX = errX
            config.plugins[self.algName].errY = errY
            task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
            exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=0)
            task.run(cat, exposure)
            return cat[0]

        # Errors are real numbers: flags should not be set.
        source = runMeasurement(1.0, 1.0)
        self.assertFalse(source.get("test_Centroider_flag"))
        self.assertFalse(source.get("test_Centroider_flag_badError"))

        # Error on X is NaN: flags should be set.
        source = runMeasurement(float('nan'), 1.0)
        self.assertTrue(source.get("test_Centroider_flag"))
        self.assertTrue(source.get("test_Centroider_flag_badError"))

        # Error on Y is NaN: flags should be set.
        source = runMeasurement(1.0, float('nan'))
        self.assertTrue(source.get("test_Centroider_flag"))
        self.assertTrue(source.get("test_Centroider_flag_badError"))

        # Error on both X and Y is NaN: flags should be set.
        source = runMeasurement(float('nan'), float('nan'))
        self.assertTrue(source.get("test_Centroider_flag"))
        self.assertTrue(source.get("test_Centroider_flag_badError"))

    def testCentroidDistance(self):
        """Test that a slight centroid movement triggers the distance error.
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeSingleFrameMeasurementConfig(plugin=self.algName)
        config.plugins[self.algName].moveX = -2
        config.plugins[self.algName].dist = 1
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=1)
        source = cat[0]
        task.run(cat, exposure)
        self.assertTrue(source.get("test_Centroider_flag"))
        self.assertTrue(source.get("test_Centroider_flag_resetToPeak"))
        self.assertEqual(source.getFootprint().getPeaks()[0].getFx(), source.get("test_Centroider_x"))

    def testCentroidOutsideFootprint(self):
        """A large centroid movement should trigger a move back to first peak.
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeSingleFrameMeasurementConfig(plugin=self.algName)
        config.plugins[self.algName].moveX = -30
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=2)
        source = cat[0]
        task.run(cat, exposure)
        self.assertTrue(source.get("test_Centroider_flag"))
        self.assertTrue(source.get("test_Centroider_flag_resetToPeak"))
        self.assertEqual(source.getFootprint().getPeaks()[0].getFx(), source.get("test_Centroider_x"))

    def testSdssCentroid(self):
        """Test the `SdssCentroid` works with the ``maxDistance`` check.
        """
        schema = self.dataset.makeMinimalSchema()
        self.algName = "base_SdssCentroid"
        config = self.makeSingleFrameMeasurementConfig(plugin=self.algName)
        config.plugins[self.algName].maxDistToPeak = .0001
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=4)
        source = cat[0]
        task.run(cat, exposure)
        self.assertTrue(source.get("base_SdssCentroid_flag"))
        self.assertTrue(source.get("base_SdssCentroid_flag_resetToPeak"))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
