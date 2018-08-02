#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
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

import numpy as np

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

    moveX = lsst.pex.config.Field(dtype=int, default=0, optional=False,
                                  doc="amount to re-position in X")
    moveY = lsst.pex.config.Field(dtype=int, default=0, optional=False,
                                  doc="amount to re-position in Y")
    dist = lsst.pex.config.Field(dtype=int, default=None, optional=False,
                                 doc="distance to allow centroid to be off")


@register("test_Centroider")
class Centroider(SingleFramePlugin):
    '''
    This is a sample Python plugin.  The flag handler for this plugin is created
    during construction, and is called using the method fail().  All plugins are
    required to implement this method, which is used to set the flags in the
    output source record if an error occurs.
    '''
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
        self.xKey = schema.addField(schema.join(name, "x"), type=np.float64)
        self.yKey = schema.addField(schema.join(name, "y"), type=np.float64)
        if self.config.dist is None:
            self.centroidChecker = lsst.meas.base.CentroidChecker(schema, name)
        else:
            self.centroidChecker = lsst.meas.base.CentroidChecker(schema, name, True, self.config.dist)

    def measure(self, measRecord, exposure):
        """
        This measure routine moves the centroid by design to create an error.
        """
        measRecord.set(self.xKey, measRecord.getX() + self.config.moveX)
        measRecord.set(self.yKey, measRecord.getY() + self.config.moveY)
        self.centroidChecker(measRecord)

    def fail(self, measRecord, error=None):
        """
        This routine responds to the standard failure call in baseMeasurement
        If the exception is a MeasurementError, the error will be passed to the
        fail method by the MeasurementFramework.
        """
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)


class FlagHandlerTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    #   Setup a configuration and datasource to be used by the plugin tests
    def setUp(self):
        self.algName = "test_Centroider"
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(100, 100), invert=False)
        self.dataset = lsst.meas.base.tests.TestDataset(bbox)
        self.dataset.addSource(flux=1E5, centroid=lsst.geom.Point2D(25, 26))

    def makeConfig(self, algName=None):
        if algName is None:
            algName = self.algName
        config = lsst.meas.base.SingleFrameMeasurementConfig()
        config.plugins = [algName]
        config.slots.centroid = None
        config.slots.apFlux = None
        config.slots.calibFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.slots.psfFlux = None
        config.slots.shape = None
        config.slots.psfShape = None
        return config

    def tearDown(self):
        del self.dataset

    def testNoError(self):
        """
        Test to be sure that the resetToPeak flag is not set when no error
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeConfig()
        config.slots.centroid = "truth"
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=0)
        task.run(cat, exposure)
        source = cat[0]
        self.assertFalse(source.get("test_Centroider_flag"))
        self.assertFalse(source.get("test_Centroider_flag_resetToPeak"))

    def testCentroidDistance(self):
        """
        test if a slight move of the centroid (but inside the footprint) will trigger
        the distance error
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeConfig()
        config.slots.centroid = "truth"
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
        """
        test if moving the centroid completely outside of the footprint will trigger
        the move back to the first peak
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeConfig()
        config.slots.centroid = "truth"
        config.plugins[self.algName].moveX = -30
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=2)
        source = cat[0]
        task.run(cat, exposure)
        self.assertTrue(source.get("test_Centroider_flag"))
        self.assertTrue(source.get("test_Centroider_flag_resetToPeak"))
        self.assertEqual(source.getFootprint().getPeaks()[0].getFx(), source.get("test_Centroider_x"))

    def testNaiveCentroid(self):
        """
        Test to be sure that NaiveCentroid centroid checker works with maxDistance check
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeConfig("base_NaiveCentroid")
        config.plugins["base_NaiveCentroid"].maxDistToPeak = .0001
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=3)
        source = cat[0]
        task.run(cat, exposure)
        self.assertTrue(source.get("base_NaiveCentroid_flag"))
        self.assertTrue(source.get("base_NaiveCentroid_flag_resetToPeak"))

    def testSdssCentroid(self):
        """
        Test to be sure that SdssCentroid centroid checker works with maxDistance check
        """
        schema = self.dataset.makeMinimalSchema()
        config = self.makeConfig("base_SdssCentroid")
        config.plugins["base_SdssCentroid"].maxDistToPeak = .0001
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
