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

from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.afw.detection import InvalidPsfError
from lsst.meas.base.tests import AlgorithmTestCase
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import FlagDefinitionList, FlagHandler


class InvalidPsfPluginConfig(SingleFramePluginConfig):
    pass


@register("test_InvalidPsfPlugin")
class InvalidPsfPlugin(SingleFramePlugin):
    ConfigClass = InvalidPsfPluginConfig

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        flagDefs = FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag()
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def measure(self, measRecord, exposure):
        raise InvalidPsfError("This record has an invalid PSF!")

    def fail(self, measRecord, error=None):
        self.flagHandler.handleFailure(measRecord)


class InvalidPsfFlagTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    def setUp(self):
        self.algName = "test_InvalidPsfPlugin"
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

    def testInvalidPsfFlag(self):
        schema = self.dataset.makeMinimalSchema()
        task = lsst.meas.base.SingleFrameMeasurementTask(schema=schema, config=self.config)
        exposure, cat = self.dataset.realize(noise=100.0, schema=schema, randomSeed=0)
        task.run(cat, exposure)

        # Check the plugin failure flag.
        self.assertTrue(cat[0][f"{self.algName}_flag"])
        # And the general invalid psf flag.
        self.assertTrue(cat[0]["base_InvalidPsf_flag"])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
