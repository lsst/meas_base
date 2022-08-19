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

import lsst.utils.tests
import lsst.geom
import lsst.meas.base
import lsst.meas.base.tests
import lsst.afw.table
from lsst.meas.base import FlagDefinitionList, FlagHandler, MeasurementError
from lsst.meas.base.tests import AlgorithmTestCase

import lsst.pex.exceptions
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin


class PythonPluginConfig(SingleFramePluginConfig):
    """Configuration for a sample plugin with a `FlagHandler`.
    """

    edgeLimit = lsst.pex.config.Field(dtype=int, default=0,
                                      doc="How close to the edge can the object be?")
    size = lsst.pex.config.Field(dtype=int, default=1,
                                 doc="size of aperture to measure around the center?")
    flux0 = lsst.pex.config.Field(dtype=float, default=None, optional=True,
                                  doc="Flux for zero mag, used to set mag if defined")


@register("test_PythonPlugin")
class PythonPlugin(SingleFramePlugin):
    """Example Python measurement plugin using a `FlagHandler`.

    This is a sample Python plugin which shows how to create and use a
    `FlagHandler`. The `FlagHandler` defines the known failures which can
    occur when the plugin is called, and should be tested after `measure` to
    detect any potential problems.

    This plugin is a very simple flux measurement algorithm which sums the
    pixel values in a square box of dimension `PythonPluginConfig.size` around
    the center point.

    Note that to properly set the error flags when a `MeasurementError` occurs,
    the plugin must implement the `fail` method as shown below. The `fail`
    method should set both the general error flag, and any specific flag as
    designated in the `MeasurementError`.

    This example also demonstrates the use of the `SafeCentroidExtractor`. The
    `SafeCentroidEextractor` and `SafeShapeExtractor` can be used to get some
    reasonable estimate of the centroid or shape in cases where the centroid
    or shape slot has failed on a particular source.
    """

    ConfigClass = PythonPluginConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        flagDefs = FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag()
        self.CONTAINS_NAN = flagDefs.add("flag_containsNan", "Measurement area contains a nan")
        self.EDGE = flagDefs.add("flag_edge", "Measurement area over edge")
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)
        self.centroidExtractor = lsst.meas.base.SafeCentroidExtractor(schema, name)
        self.instFluxKey = schema.addField(name + "_instFlux", "F", doc="flux")
        self.magKey = schema.addField(name + "_mag", "F", doc="mag")

    def measure(self, measRecord, exposure):
        """Perform measurement.

        The measure method is called by the measurement framework when
        `run` is called. If a `MeasurementError` is raised during this method,
        the `fail` method will be called to set the error flags.
        """
        # Call the SafeCentroidExtractor to get a centroid, even if one has
        # not been supplied by the centroid slot. Normally, the centroid is
        # supplied by the centroid slot, but if that fails, the footprint is
        # used as a fallback.  If the fallback is needed, the fail flag will
        # be set on this record.
        center = self.centroidExtractor(measRecord, self.flagHandler)

        # create a square bounding box of size = config.size around the center
        centerPoint = lsst.geom.Point2I(int(center.getX()), int(center.getY()))
        bbox = lsst.geom.Box2I(centerPoint, lsst.geom.Extent2I(1, 1))
        bbox.grow(self.config.size)

        # If the measurement box falls outside the exposure, raise the edge
        # MeasurementError
        if not exposure.getBBox().contains(bbox):
            raise MeasurementError(self.EDGE.doc, self.EDGE.number)

        # Sum the pixels inside the bounding box
        instFlux = lsst.afw.image.ImageF(exposure.getMaskedImage().getImage(), bbox).getArray().sum()
        measRecord.set(self.instFluxKey, instFlux)

        # If there was a NaN inside the bounding box, the instFlux will still
        # be NaN
        if np.isnan(instFlux):
            raise MeasurementError(self.CONTAINS_NAN.doc, self.CONTAINS_NAN.number)

        if self.config.flux0 is not None:
            if self.config.flux0 == 0:
                raise ZeroDivisionError("self.config.flux0 is zero in divisor")
            mag = -2.5 * np.log10(instFlux/self.config.flux0)
            measRecord.set(self.magKey, mag)

    def fail(self, measRecord, error=None):
        """Handle measurement failures.

        If the exception is a `MeasurementError`, the error will be passed to
        the fail method by the measurement Framework. If ``error`` is not
        `None`, ``error.cpp`` should correspond to a specific error and the
        appropriate error flag will be set.
        """
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)


class FlagHandlerTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    # Setup a configuration and datasource to be used by the plugin tests

    def setUp(self):
        self.algName = "test_PythonPlugin"
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(bbox)
        self.dataset.addSource(instFlux=1E5, centroid=lsst.geom.Point2D(25, 26))
        config = lsst.meas.base.SingleFrameMeasurementConfig()
        config.plugins = [self.algName]
        config.slots.centroid = None
        config.slots.apFlux = None
        config.slots.calibFlux = None
        config.slots.gaussianFlux = None
        config.slots.modelFlux = None
        config.slots.psfFlux = None
        config.slots.shape = None
        config.slots.psfShape = None
        self.config = config

    def tearDown(self):
        del self.config
        del self.dataset

    def testFlagHandler(self):
        """Test creation and invocation of `FlagHander`.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()

        # This is a FlagDefinition structure like a plugin might have
        flagDefs = FlagDefinitionList()
        FAILURE = flagDefs.addFailureFlag()
        FIRST = flagDefs.add("1st error", "this is the first failure type")
        SECOND = flagDefs.add("2nd error", "this is the second failure type")
        fh = FlagHandler.addFields(schema, "test", flagDefs)
        # Check to be sure that the FlagHandler was correctly initialized
        for index in range(len(flagDefs)):
            self.assertEqual(flagDefs.getDefinition(index).name, fh.getFlagName(index))

        catalog = lsst.afw.table.SourceCatalog(schema)

        # Now check to be sure that all of the known failures set the bits
        # correctly
        record = catalog.addNew()
        fh.handleFailure(record)
        self.assertTrue(fh.getValue(record, FAILURE.number))
        self.assertFalse(fh.getValue(record, FIRST.number))
        self.assertFalse(fh.getValue(record, SECOND.number))
        record = catalog.addNew()

        error = MeasurementError(FAILURE.doc, FAILURE.number)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FAILURE.number))
        self.assertFalse(fh.getValue(record, FIRST.number))
        self.assertFalse(fh.getValue(record, SECOND.number))

        record = catalog.addNew()
        error = MeasurementError(FIRST.doc, FIRST.number)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FAILURE.number))
        self.assertTrue(fh.getValue(record, FIRST.number))
        self.assertFalse(fh.getValue(record, SECOND.number))

        record = catalog.addNew()
        error = MeasurementError(SECOND.doc, SECOND.number)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FAILURE.number))
        self.assertFalse(fh.getValue(record, FIRST.number))
        self.assertTrue(fh.getValue(record, SECOND.number))

    def testNoFailureFlag(self):
        """Test with no failure flag.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()

        # This is a FlagDefinition structure like a plugin might have
        flagDefs = FlagDefinitionList()
        FIRST = flagDefs.add("1st error", "this is the first failure type")
        SECOND = flagDefs.add("2nd error", "this is the second failure type")
        fh = FlagHandler.addFields(schema, "test", flagDefs)
        # Check to be sure that the FlagHandler was correctly initialized
        for index in range(len(flagDefs)):
            self.assertEqual(flagDefs.getDefinition(index).name, fh.getFlagName(index))

        catalog = lsst.afw.table.SourceCatalog(schema)

        # Now check to be sure that all of the known failures set the bits
        # correctly
        record = catalog.addNew()
        fh.handleFailure(record)
        self.assertFalse(fh.getValue(record, FIRST.number))
        self.assertFalse(fh.getValue(record, SECOND.number))
        record = catalog.addNew()

        record = catalog.addNew()
        error = MeasurementError(FIRST.doc, FIRST.number)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, FIRST.number))
        self.assertFalse(fh.getValue(record, SECOND.number))

        record = catalog.addNew()
        error = MeasurementError(SECOND.doc, SECOND.number)
        fh.handleFailure(record, error.cpp)
        self.assertFalse(fh.getValue(record, FIRST.number))
        self.assertTrue(fh.getValue(record, SECOND.number))

    # This and the following tests using the toy plugin, and demonstrate how
    # the FlagHandler is used.

    def testPluginNoError(self):
        """Test that the sample plugin can be run without errors.
        """
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=0)
        task.run(cat, exposure)
        source = cat[0]
        self.assertFalse(source.get(self.algName + "_flag"))
        self.assertFalse(source.get(self.algName + "_flag_containsNan"))
        self.assertFalse(source.get(self.algName + "_flag_edge"))

    def testPluginUnexpectedError(self):
        """Test that unexpected non-fatal errors set the failure flag.

        An unexpected error is a non-fatal error which is not caught by the
        algorithm itself.  However, such errors are caught by the measurement
        framework in task.run, and result in the failure flag being set, but
        no other specific flags
        """
        self.config.plugins[self.algName].flux0 = 0.0  # divide by zero
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=1)
        task.log.setLevel(task.log.FATAL)
        task.run(cat, exposure)
        source = cat[0]
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertFalse(source.get(self.algName + "_flag_containsNan"))
        self.assertFalse(source.get(self.algName + "_flag_edge"))

    def testPluginContainsNan(self):
        """Test that the ``containsNan`` error can be triggered.
        """
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=2)
        source = cat[0]
        exposure.getMaskedImage().getImage().getArray()[int(source.getY()), int(source.getX())] = np.nan
        task.run(cat, exposure)
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertTrue(source.get(self.algName + "_flag_containsNan"))
        self.assertFalse(source.get(self.algName + "_flag_edge"))

    def testPluginEdgeError(self):
        """Test that the ``edge`` error can be triggered.
        """
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=3)
        # Set the size large enough to trigger the edge error
        self.config.plugins[self.algName].size = exposure.getDimensions()[1]//2
        task.log.setLevel(task.log.FATAL)
        task.run(cat, exposure)
        source = cat[0]
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertFalse(source.get(self.algName + "_flag_containsNan"))
        self.assertTrue(source.get(self.algName + "_flag_edge"))

    def testSafeCentroider(self):
        """Test `SafeCentroidExtractor` correctly runs and sets flags.
        """
        # Normal case should use the centroid slot to get the center, which
        # should succeed
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        task.log.setLevel(task.log.FATAL)
        exposure, cat = self.dataset.realize(noise=0.0, schema=schema, randomSeed=4)
        source = cat[0]
        task.run(cat, exposure)
        self.assertFalse(source.get(self.algName + "_flag"))
        instFlux = source.get("test_PythonPlugin_instFlux")
        self.assertFalse(np.isnan(instFlux))

        # If one of the center coordinates is nan and the centroid slot error
        # flag has not been set, the SafeCentroidExtractor will fail.
        source.set('truth_x', np.nan)
        source.set('truth_flag', False)
        source.set("test_PythonPlugin_instFlux", np.nan)
        source.set(self.algName + "_flag", False)
        task.run(cat, exposure)
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertTrue(np.isnan(source.get("test_PythonPlugin_instFlux")))

        # But if the same conditions occur and the centroid slot error flag is
        # set to true, the SafeCentroidExtractor will succeed and the
        # algorithm will complete.  However, the failure flag will also be
        # set.
        source.set('truth_x', np.nan)
        source.set('truth_flag', True)
        source.set("test_PythonPlugin_instFlux", np.nan)
        source.set(self.algName + "_flag", False)
        task.run(cat, exposure)
        self.assertTrue(source.get(self.algName + "_flag"))
        self.assertEqual(source.get("test_PythonPlugin_instFlux"), instFlux)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
