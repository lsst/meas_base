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
""" Unit tests for Python Plugin FlagHandlers and Sample Plugin Example."""
import numpy as np
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
from lsst.meas.base.flagDecorator import addFlagHandler


class PythonPluginConfig(SingleFramePluginConfig):
    """
    Configuration for Sample Plugin with a FlagHandler.
    """

    edgeLimit = lsst.pex.config.Field(dtype=int, default=0, optional=False,
                                      doc="How close to the edge can the object be?")
    size = lsst.pex.config.Field(dtype=int, default=1, optional=False,
                                 doc="size of aperture to measure around the center?")
    flux0 = lsst.pex.config.Field(dtype=float, default=None, optional=False,
                                  doc="Flux for zero mag, used to set mag if defined")


@register("test_PythonPlugin")
@addFlagHandler(("flag", "General Failure error"),
                ("flag_containsNan", "Measurement area contains a nan"),
                ("flag_edge", "Measurement area over edge"))
class PythonPlugin(SingleFramePlugin):
    """
    This is a sample Python plugin which shows how to create and use a FlagHandler.
    The FlagHandler defines the known failures which can occur when the plugin is
    called, and should be tested after measure() to detect any potential problems.

    This plugin is a very simple flux measurement algorithm which sums
    the pixel values in a square box of dimension config.size around the center point.

    The flagHandler for this plugin is created during construction as a result of the
    addFlagHandler decorator, which adds a general failure flag and two specific flags
    (containsNan and edge) to the schema. The class variable ErrEnum is also added.
    This is an enumeration of the 3 flags.

    Note that to properly set the error flags when a MeasurementError occurs, the plugin
    must implement the fail() method as shown below. The fail method should set both
    the general error flag, and any specific flag (such as ErrEnum.flag_edge or
    ErrEnum.flag_containsNan) as designated in the MeasurementError.

    This example also demonstrates the use of the SafeCentroidExtractor. The
    SafeCentroidEextractor and SafeShapeExtractor can be used to get some reasonable
    estimate of the centroid or shape in cases where the centroid or shape slot
    has failed on a particular source.
    """

    ConfigClass = PythonPluginConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    # Note that because of the decorator, a flagHander and ErrEnum are also created
    # and the flags are added to the schema.
    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.centroidExtractor = lsst.meas.base.SafeCentroidExtractor(schema, name)
        self.fluxKey = schema.addField(name + "_flux", "F", doc="flux")
        self.magKey = schema.addField(name + "_mag", "F", doc="mag")

    def measure(self, measRecord, exposure):
        """
        The measure method is called by the measurement framework when task.run is called
        If a MeasurementError is raised during this method, the fail() method will be
        called to set the error flags.
        """
        # Call the SafeCentroidExtractor to get a centroid, even if one has
        # not been supplied by the centroid slot. Normally, the centroid is supplied
        # by the centroid slot, but if that fails, the footprint is used as a fallback.
        # If the fallback is needed, the fail flag will be set on this record.
        center = self.centroidExtractor(measRecord, self.flagHandler)

        # create a square bounding box of size = config.size around the center
        centerPoint = lsst.afw.geom.Point2I(int(center.getX()), int(center.getY()))
        bbox = lsst.afw.geom.Box2I(centerPoint, lsst.afw.geom.Extent2I(1, 1))
        bbox.grow(self.config.size)

        # If the measurement box falls outside the exposure, raise the edge MeasurementError
        if not exposure.getBBox().contains(bbox):
            raise MeasurementError(self.flagHandler.getDefinition(self.ErrEnum.flag_edge).doc,
                                   PythonPlugin.ErrEnum.flag_edge)

        # Sum the pixels inside the bounding box
        flux = lsst.afw.image.ImageF(exposure.getMaskedImage().getImage(), bbox).getArray().sum()
        measRecord.set(self.fluxKey, flux)

        # If there was a nan inside the bounding box, the flux will still be nan
        if np.isnan(flux):
            raise MeasurementError(self.flagHandler.getDefinition(self.ErrEnum.flag_containsNan).doc,
                                   PythonPlugin.ErrEnum.flag_containsNan)

        if self.config.flux0 is not None:
            if self.config.flux0 == 0:
                raise ZeroDivisionError("self.config.flux0 is zero in divisor")
            mag = -2.5 * np.log10(flux/self.config.flux0)
            measRecord.set(self.magKey, mag)

    def fail(self, measRecord, error=None):
        """
        This routine responds to the standard failure call in baseMeasurement
        If the exception is a MeasurementError, the error will be passed to the
        fail method by the MeasurementFramework. If error is not none, error.cpp
        should correspond to a specific error in the ErrEnum, and the appropriate
        error flag will be set.
        """
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)


class FlagHandlerTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    # Setup a configuration and datasource to be used by the plugin tests
    def setUp(self):
        self.algName = "test_PythonPlugin"
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0), lsst.afw.geom.Point2I(100, 100))
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

    def testFlagHandler(self):
        """
        Standalone test to create a flaghandler and call it
        This is not a real world example, just a simple unit test
        """
        #control = lsst.meas.base.GaussianCentroidControl()
        #alg = lsst.meas.base.GaussianCentroidAlgorithm
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        #plugin = alg(control, 'test', schema)
        #cat = lsst.afw.table.SourceCatalog(schema)
        #subSchema = schema["test"]

        # This is a FlagDefinition structure like a plugin might have
        FAILURE = 0
        FIRST = 1
        SECOND = 2
        flagDefs = [FlagDefinition("General Failure", "general failure error"),
                    FlagDefinition("1st error", "this is the first failure type"),
                    FlagDefinition("2nd error", "this is the second failure type")
                    ]
        fh = FlagHandler.addFields(schema, "test", FlagDefinitionVector(flagDefs))

        # Check to be sure that the FlagHandler was correctly initialized
        for index, flagDef in enumerate(flagDefs):
            self.assertEqual(flagDef.name, fh.getDefinition(index).name)
            self.assertEqual(flagDef.doc, fh.getDefinition(index).doc)

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

    # This and the following tests using the toy plugin, and demonstrate how
    # the flagHandler is used.

    def testPluginNoError(self):
        """
        Test that the sample plugin can be run without errors
        """
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        task.run(cat, exposure)
        source = cat[0]
        self.assertFalse(source.get(self.algName + "_flag"))
        self.assertFalse(source.get(self.algName + "_flag_containsNan"))
        self.assertFalse(source.get(self.algName + "_flag_edge"))

    def testPluginUnexpectedError(self):
        """
        An unexpected error is a non-fatal error which is not caught by the algorithm itself.
        However, such errors are caught by the measurement framework in task.run, and result
        in the failure flag being set, but no other specific flags
        """
        self.config.plugins[self.algName].flux0 = 0.0     # this causes a divide by zero
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        task.log.setThreshold(task.log.FATAL)
        task.run(cat, exposure)
        source = cat[0]
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertFalse(source.get(self.algName + "_flag_containsNan"))
        self.assertFalse(source.get(self.algName + "_flag_edge"))

    def testPluginContainsNan(self):
        """
        Test that the containsNan error can be triggered.
        """
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        source = cat[0]
        exposure.getMaskedImage().getImage().getArray()[int(source.getY()), int(source.getX())] = np.nan
        task.run(cat, exposure)
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertTrue(source.get(self.algName + "_flag_containsNan"))
        self.assertFalse(source.get(self.algName + "_flag_edge"))

    def testPluginEdgeError(self):
        """
        Test that the Edge error can be triggered.
        """
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema)
        # Set the size large enough to trigger the edge error
        self.config.plugins[self.algName].size = exposure.getDimensions()[1]//2
        task.log.setThreshold(task.log.FATAL)
        task.run(cat, exposure)
        source = cat[0]
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertFalse(source.get(self.algName + "_flag_containsNan"))
        self.assertTrue(source.get(self.algName + "_flag_edge"))

    def testSafeCentroider(self):
        """
        Test that the SafeCentroidExtractor works and sets the flags correctly
        """
        # Normal case should use the centroid slot to get the center, which should succeed
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        task.log.setThreshold(task.log.FATAL)
        exposure, cat = self.dataset.realize(noise=0.0, schema=schema)
        source = cat[0]
        task.run(cat, exposure)
        self.assertFalse(source.get(self.algName + "_flag"))
        flux = source.get("test_PythonPlugin_flux")
        self.assertFalse(np.isnan(flux))

        # If one of the center coordinates is nan and the centroid slot error flag has
        # not been set, the SafeCentroidExtractor will fail.
        source.set('truth_x', np.nan)
        source.set('truth_flag', False)
        source.set("test_PythonPlugin_flux", np.nan)
        source.set(self.algName + "_flag", False)
        task.run(cat, exposure)
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertTrue(np.isnan(source.get("test_PythonPlugin_flux")))

        # But if the same conditions occur and the centroid slot error flag is set
        # to true, the SafeCentroidExtractor will succeed and the algorithm will complete.
        # However, the failure flag will also be set.
        source.set('truth_x', np.nan)
        source.set('truth_flag', True)
        source.set("test_PythonPlugin_flux", np.nan)
        source.set(self.algName + "_flag", False)
        task.run(cat, exposure)
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertEqual(source.get("test_PythonPlugin_flux"), flux)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
