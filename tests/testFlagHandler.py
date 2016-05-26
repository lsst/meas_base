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
import os
import unittest

import lsst.utils.tests
import lsst.meas.base
import lsst.meas.base.tests
import lsst.afw.table
from lsst.meas.base.baseLib import MeasurementError
from lsst.meas.base import FlagDefinition, FlagDefinitionVector, FlagHandler
from lsst.meas.base.tests import (AlgorithmTestCase)

import lsst.pex.exceptions
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base.baseLib import MeasurementError
from lsst.meas.base import FlagDefinition, FlagDefinitionVector, FlagHandler
from lsst.meas.base.flagDecorator import addFlagHandler


class PythonPluginConfig(SingleFramePluginConfig):

    failureType = lsst.pex.config.Field(dtype=int, default=None, optional=False,
                                  doc="A failure mode to test")

@register("test_PythonPlugin")
@addFlagHandler(("flag", "General Failure error"), ("flag_error1","First type of Failure occured."),
                ("flag_error2", "Second type of failure occured."))
class PythonPlugin(SingleFramePlugin):
    '''
    This is a sample Python plugin.  The flag handler for this plugin is created
    during construction, and is called using the method fail().  All plugins are
    required to implement this method, which is used to set the flags in the
    output source record if an error occurs.
    '''
    ConfigClass = PythonPluginConfig
    # Class variables ErrEnum and FLAGDEFS are added by the decorator

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        # The instance variable flagHandler is added by the decorator

    #   This is a measure routine which does nothing except to raise Exceptions
    #   as requested by the caller. Errors normally don't occur unless there is
    #   something wrong in the inputs, or if there is an error during the measurement
    def measure(self, measRecord, exposure):
        if not self.config.failureType is None:
            if self.config.failureType == PythonPlugin.ErrEnum.flag_error1:
                raise MeasurementError(self.flagHandler.getDefinition(PythonPlugin.ErrEnum.flag_error1).doc,
                    PythonPlugin.ErrEnum.flag_error1)
            if self.config.failureType == PythonPlugin.ErrEnum.flag_error2:
                raise MeasurementError(self.flagHandler.getDefinition(PythonPlugin.ErrEnum.flag_error2).doc,
                    PythonPlugin.ErrEnum.flag_error2)
            raise RuntimeError("An unexpected error occurred")

    #   This routine responds to the standard failure call in baseMeasurement
    #   If the exception is a MeasurementError, the error will be passed to the
    #   fail method by the MeasurementFramework.
    def fail(self, measRecord, error=None):
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)

class FlagHandlerTestCase(AlgorithmTestCase):

    #   Setup a configuration and datasource to be used by the plugin tests
    def setUp(self):
        self.algName = "test_PythonPlugin"
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0,0), lsst.afw.geom.Point2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(bbox)
        self.dataset.addSource(flux=1E5, centroid=lsst.afw.geom.Point2D(25, 26))
        config = lsst.meas.base.SingleFrameMeasurementConfig()
        config.plugins = [self.algName]
        config.slots.centroid = None
        config.slots.apFlux = None
        config.slots.calibFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.slots.psfFlux = None
        config.slots.shape = None
        self.config = config

    def tearDown(self):
        del self.config
        del self.dataset

    #   Standalone test to create a flaghandler and call it
    #   This is not a real world example, just a simple unit test
    def testFlagHandler(self):
        control = lsst.meas.base.GaussianCentroidControl()
        alg = lsst.meas.base.GaussianCentroidAlgorithm
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        plugin = alg(control, 'test', schema)
        cat = lsst.afw.table.SourceCatalog(schema)
        subSchema = schema["test"]

        # This is a FlagDefinition structure like a plugin might have
        FAILURE = 0
        FIRST = 1
        SECOND = 2
        flagDefs = [ FlagDefinition("General Failure", "general failure error"),
            FlagDefinition("1st error", "this is the first failure type"),
            FlagDefinition("2nd error", "this is the second failure type")
        ]
        fh = FlagHandler.addFields(schema, "test",
            FlagDefinitionVector(flagDefs))

        # Check to be sure that the FlagHandler was correctly initialized
        for index, flagDef in enumerate(flagDefs):
           assert(flagDef.name == fh.getDefinition(index).name)
           assert(flagDef.doc == fh.getDefinition(index).doc)

        catalog = lsst.afw.table.SourceCatalog(schema)

        # Now check to be sure that all of the known failures set the bits correctly
        record = catalog.addNew()
        fh.handleFailure(record)
        self.assertTrue(fh.getValue(record, FAILURE))
        self.assertFalse(fh.getValue(record, FIRST))
        self.assertFalse(fh.getValue(record, SECOND))
        record = catalog.addNew()

        error = MeasurementError(fh.getDefinition(FAILURE).doc, FAILURE)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FAILURE))
        self.assertFalse(fh.getValue(record, FIRST))
        self.assertFalse(fh.getValue(record, SECOND))

        record = catalog.addNew()
        error = MeasurementError(fh.getDefinition(FIRST).doc, FIRST)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FAILURE))
        self.assertTrue(fh.getValue(record, FIRST))
        self.assertFalse(fh.getValue(record, SECOND))

        record = catalog.addNew()
        error = MeasurementError(fh.getDefinition(SECOND).doc, SECOND)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FAILURE))
        self.assertFalse(fh.getValue(record, FIRST))
        self.assertTrue(fh.getValue(record, SECOND))

    def testNoError(self):
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        task.run(cat, exposure)
        source = cat[0]
        self.assertEqual(source.get(self.algName + "_flag"), False)
        self.assertEqual(source.get(self.algName + "_flag_error1"), False)
        self.assertEqual(source.get(self.algName + "_flag_error2"), False)

    def testUnexpectedError(self):
        self.config.plugins[self.algName].failureType = -1     # any unknown error type will do
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        task.log.setThreshold(task.log.FATAL)
        task.run(cat, exposure)
        source = cat[0]
        self.assertEqual(source.get(self.algName + "_flag"), True)
        self.assertEqual(source.get(self.algName + "_flag_error1"), False)
        self.assertEqual(source.get(self.algName + "_flag_error2"), False)

    def testError1(self):
        self.config.plugins[self.algName].failureType = PythonPlugin.ErrEnum.flag_error1
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        task.run(cat, exposure)
        source = cat[0]
        self.assertEqual(source.get(self.algName + "_flag"), True)
        self.assertEqual(source.get(self.algName + "_flag_error1"), True)
        self.assertEqual(source.get(self.algName + "_flag_error2"), False)

    def testError2(self):
        self.config.plugins[self.algName].failureType = PythonPlugin.ErrEnum.flag_error2
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=0.0, schema=schema)
        task.run(cat, exposure)
        source = cat[0]
        self.assertEqual(source.get(self.algName + "_flag"), True)
        self.assertEqual(source.get(self.algName + "_flag_error1"), False)
        self.assertEqual(source.get(self.algName + "_flag_error2"), True)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(FlagHandlerTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
