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
import os
import numpy

import lsst.geom
import lsst.afw.table
import lsst.daf.base
import lsst.meas.base
import lsst.utils.tests
from lsst.meas.base.tests import (AlgorithmTestCase, )
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base.forcedMeasurement import ForcedPlugin
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import FlagDefinitionList, FlagHandler, MeasurementError

ROOT = os.path.abspath(os.path.dirname(__file__))


class LoggingPluginConfig(SingleFramePluginConfig):
    """Configuration for sample plugin.
    """
    pass


@register("test_LoggingPlugin")
class LoggingPlugin(SingleFramePlugin):
    """Sample Python plugin which has an associated log name.

    Notes
    -----
    The log name is provided to the plugin by the measurement task which is
    running it. This requires that the `hasLogName` attribute must be a member
    of the plugin class, and it must be `True`.
    """
    hasLogName = True
    ConfigClass = LoggingPluginConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    # The initializer for the class must accept an optional logName parameter.
    def __init__(self, config, name, schema, metadata, logName=None):
        SingleFramePlugin.__init__(self, config, name, schema, metadata, logName=logName)
        flagDefs = FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag()
        self.CONTAINS_NAN = flagDefs.add("flag_containsNan", "Measurement area contains a nan")
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)
        self.instFluxKey = schema.addField(name + "_instFlux", "F", doc="flux")

    def measure(self, measRecord, exposure):
        """Perform measurement.

        Notes
        -----
        The `measure` method is called by the measurement framework when `run`
        is called. If a `MeasurementError` is raised during this method, the
        `fail` method will be called to set the error flags.
        """
        lsst.log.Log.getLogger(self.getLogName()).info("%s plugin measuring."%(self.name,))
        # Sum the pixels inside the bounding box
        centerPoint = lsst.geom.Point2I(int(measRecord.getX()), int(measRecord.getY()))
        bbox = lsst.geom.Box2I(centerPoint, lsst.geom.Extent2I(1, 1))
        instFlux = lsst.afw.image.ImageF(exposure.getMaskedImage().getImage(), bbox).getArray().sum()
        measRecord.set(self.instFluxKey, instFlux)

        # If there was a NaN inside the bounding box, the instFlux will still
        # be NaN
        if numpy.isnan(instFlux):
            raise MeasurementError(self.CONTAINS_NAN.doc, self.CONTAINS_NAN.number)

    def fail(self, measRecord, error=None):
        """Handle measurement failures.

        Notes
        -----
        If measurement raises a `MeasurementError`, the error will be passed
        to the fail method by the measurement framework. If the error is not
        `None`, ``error.cpp`` should correspond to a specific error and the
        appropriate error flag will be set.
        """
        if error is None:
            self.flagHandler.handleFailure(measRecord)
        else:
            self.flagHandler.handleFailure(measRecord, error.cpp)


def directLog(log, file=None):
    """Direct the log given to a file or to the console if ``file`` is `None`.
    """
    props = "log4j.rootLogger=INFO, FA\n"
    if file is None:
        props += "log4j.appender.FA=ConsoleAppender\n"
    else:
        props += "log4j.appender.FA=FileAppender\n"
        props += "log4j.appender.FA.Append=false\n"
        props += "log4j.appender.FA.file=%s\n"%(file,)
        props += "log4j.appender.FA.Append=false\n"
    props += "log4j.appender.FA.layout=PatternLayout\n"
    props += "log4j.appender.FA.layout.ConversionPattern=%d{yyyy-MM-dd HH:mm:ss.SSS} %p %c %m %X%n\n"
    props += "log4j.logger.main.a=DEBUG\n"
    log.configure_prop(props)


class RegisteredPluginsTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    """Test all registered Plugins to see if their logName is set as expected.

    Those which have the ``hasLogName=True`` attribute will have a ``logName``
    parameter passed to their ``__init__``, and should set the internal
    ``_logName`` attribute.  If they are wrapped C++ algorithms, the
    `getLogName` should also return same ``logName`` as the plugin.
    """
    def testSingleFramePlugins(self):
        center = lsst.geom.Point2D(50, 50)
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Extent2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        dataset.addSource(1000000.0, center)
        registry = SingleFramePlugin.registry
        dependencies = registry.keys()
        task = self.makeSingleFrameMeasurementTask("base_SdssCentroid", dependencies=dependencies)
        exposure, catalog = dataset.realize(noise=100.0, schema=task.schema, randomSeed=0)
        task.log.setLevel(lsst.log.ERROR)
        task.run(catalog, exposure)
        for pluginName in dependencies:
            plugin = task.plugins[pluginName]
            if hasattr(plugin, "hasLogName") and plugin.hasLogName:
                self.assertEqual(plugin.getLogName(), task.getPluginLogName(pluginName))
                # if the plugin is cpp, check the cpp Algorithm as well
                if hasattr(plugin, "cpp"):
                    self.assertEqual(plugin.cpp.getLogName(), plugin.getLogName())
            else:
                self.assertEqual(plugin.getLogName(), None)

    def testForcedPlugins(self):
        # Test all the ForcedPlugins registered to see if their logName is set
        # as expected.
        center = lsst.geom.Point2D(50, 50)
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Extent2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        dataset.addSource(1000000.0, center)
        registry = ForcedPlugin.registry
        dependencies = registry.keys()

        task = self.makeForcedMeasurementTask("base_SdssCentroid", dependencies=dependencies)
        measWcs = dataset.makePerturbedWcs(dataset.exposure.getWcs(), randomSeed=1)
        measDataset = dataset.transform(measWcs)
        exposure, truthCatalog = measDataset.realize(10.0, measDataset.makeMinimalSchema(), randomSeed=1)
        refCat = dataset.catalog
        refWcs = dataset.exposure.getWcs()
        measCat = task.generateMeasCat(exposure, refCat, refWcs)
        task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)

        task.log.setLevel(lsst.log.ERROR)
        task.run(measCat, exposure, refCat, refWcs)
        for pluginName in dependencies:
            plugin = task.plugins[pluginName]
            if hasattr(plugin, "hasLogName") and plugin.hasLogName:
                self.assertEqual(plugin.getLogName(), task.getPluginLogName(pluginName))
                # if the plugin is cpp, check the cpp Algorithm as well
                if hasattr(plugin, "cpp"):
                    self.assertEqual(plugin.cpp.getLogName(), task.log.getName() + "." + pluginName)
            else:
                self.assertEqual(plugin.getLogName(), None)


class LoggingPythonTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    """Test one C++ and one Python plugin which have hasLogName=True.
    """
    def setUp(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(bbox)
        self.dataset.addSource(instFlux=1E5, centroid=lsst.geom.Point2D(25, 25))
        config = lsst.meas.base.SingleFrameMeasurementConfig()
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

    def testLoggingPythonPlugin(self):
        algName = "test_LoggingPlugin"
        schema = self.dataset.makeMinimalSchema()
        self.config.plugins = [algName]
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        #  test that the plugin's logName has been propagated to the plugin
        self.assertTrue(task.plugins[algName].getLogName(), task.getPluginLogName(algName))
        log = lsst.log.Log.getLogger(task.getPluginLogName(algName))
        with lsst.utils.tests.getTempFilePath(".log") as pluginLogName:
            directLog(log, pluginLogName)
            exposure, cat = self.dataset.realize(noise=0.0, schema=schema, randomSeed=2)
            task.run(cat, exposure)
            directLog(log, None)
            # direct back to console, closing log files
            with open(pluginLogName) as fin:
                lines = fin.read()
        # test that the sample plugin has correctly logged to where we
        # expected it to.
        self.assertTrue(lines.find("measuring") >= 0)

    def testLoggingCppPlugin(self):
        # PsfFlux is known to log an ``ERROR`` if a Psf is not attached
        algName = "base_PsfFlux"
        self.config.plugins = [algName]

        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        log = lsst.log.Log.getLogger(task.getPluginLogName(algName))
        log.setLevel(lsst.log.ERROR)

        # test that the plugin's logName has been propagated to the plugin
        self.assertTrue(task.plugins[algName].getLogName(), task.getPluginLogName(algName))
        self.assertTrue(task.plugins[algName].cpp.getLogName(), task.getPluginLogName(algName))
        with lsst.utils.tests.getTempFilePath(".log") as pluginLogName:
            directLog(log, pluginLogName)
            exposure, cat = self.dataset.realize(noise=0.0, schema=schema, randomSeed=3)
            exposure.setPsf(None)
            # This call throws an error, so be prepared for it
            try:
                task.run(cat, exposure)
            except Exception:
                pass
            directLog(log, None)
            # direct back to console, closing log files
            with open(pluginLogName) as fin:
                lines = fin.read()
        # test that the sample plugin has correctly logged to where we
        # expected it to.
        self.assertTrue(lines.find("ERROR") >= 0)


class SingleFrameTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        # object in corner to trigger EDGE error
        self.center = lsst.geom.Point2D(5, 5)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                    lsst.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(1000000.0, self.center)
        self.task = self.makeSingleFrameMeasurementTask("base_SdssCentroid")
        self.log = lsst.log.Log.getLogger(self.task.getPluginLogName("base_SdssCentroid"))
        self.exposure, self.catalog = self.dataset.realize(10.0, self.task.schema, randomSeed=4)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset
        del self.task
        del self.log
        del self.exposure
        del self.catalog

    def testSeparatePluginLogs(self):
        """Check that the task log and the plugin log are truly separate.
        """
        taskLogName = os.path.join(ROOT, 'testSeparatePluginLogs-task.log')
        directLog(self.task.log, taskLogName)
        self.task.log.info("Testing")
        with lsst.utils.tests.getTempFilePath(".log") as pluginLogName:
            directLog(self.log, pluginLogName)
            self.log.setLevel(lsst.log.DEBUG)
            self.task.run(self.catalog, self.exposure)
            # direct back to console, closing log files
            directLog(self.log, None)
            directLog(self.task.log, None)
            with open(taskLogName) as fin:
                lines = fin.read()
            os.unlink(taskLogName)
            self.assertTrue(lines.find("Testing") >= 0)
            with open(pluginLogName) as fin:
                lines = fin.read()
        self.assertTrue(lines.find("MeasurementError") >= 0)

    def testSetPluginLevel(self):
        """Test setting the plugin log level.

        Specifically, we set it to the ``ERROR`` level.
        """
        with lsst.utils.tests.getTempFilePath(".log") as pluginLogName:
            directLog(self.log, pluginLogName)
            self.log.setLevel(lsst.log.ERROR)
            self.task.run(self.catalog, self.exposure)
            # direct back to console, closing log files
            directLog(self.log, None)
            with open(pluginLogName) as fin:
                lines = fin.read()
        self.assertTrue(lines.find("MeasurementError") < 0)


class ForcedTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        # object in corner to trigger EDGE error
        self.center = lsst.geom.Point2D(0, 0)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                    lsst.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(1000000.0, self.center)
        self.task = self.makeForcedMeasurementTask("base_SdssCentroid")
        self.log = lsst.log.Log.getLogger(self.task.getPluginLogName("base_SdssCentroid"))
        measWcs = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs(), randomSeed=5)
        measDataset = self.dataset.transform(measWcs)
        self.exposure, truthCatalog = measDataset.realize(10.0, measDataset.makeMinimalSchema(), randomSeed=5)
        self.refCat = self.dataset.catalog
        self.refWcs = self.dataset.exposure.getWcs()
        self.measCat = self.task.generateMeasCat(self.exposure, self.refCat, self.refWcs)
        self.task.attachTransformedFootprints(self.measCat, self.refCat, self.exposure, self.refWcs)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset
        del self.task
        del self.log
        del self.exposure
        del self.measCat
        del self.refCat
        del self.refWcs

    def testSeparatePluginLog(self):
        """Check that the task log and the plugin log are truly separate.
        """
        taskLogName = os.path.join(ROOT, 'testSeparatePluginLog-task.log')
        directLog(self.task.log, taskLogName)
        self.task.log.info("Testing")
        with lsst.utils.tests.getTempFilePath(".log") as pluginLogName:
            directLog(self.log, pluginLogName)
            self.log.setLevel(lsst.log.DEBUG)
            self.task.run(self.measCat, self.exposure, self.refCat, self.refWcs)
            # direct back to console, closing log files
            directLog(self.log, None)
            directLog(self.task.log, None)
            with open(taskLogName) as fin:
                lines = fin.read()
            os.unlink(taskLogName)
            self.assertTrue(lines.find("Testing") >= 0)
            with open(pluginLogName) as fin:
                lines = fin.read()
        self.assertTrue(lines.find("MeasurementError") >= 0)

    def testSetPluginLevel(self):
        """Test setting the plugin log level.

        Specifically, we set it to the ``ERROR`` level.
        """
        with lsst.utils.tests.getTempFilePath(".log") as pluginLogName:
            directLog(self.log, pluginLogName)
            self.log.setLevel(lsst.log.ERROR)
            self.task.run(self.measCat, self.exposure, self.refCat, self.refWcs)
            # direct back to console, closing log files
            directLog(self.log, None)
            with open(pluginLogName) as fin:
                lines = fin.read()
        self.assertTrue(lines.find("MeasurementError") < 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
