#
# LSST Data Management System
# Copyright 2008-2018 AURA/LSST.
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

import unittest

import numpy as np

import lsst.afw.geom
import lsst.afw.image
import lsst.utils.tests
from lsst.meas.base import LocalBackgroundAlgorithm
from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)


class LocalBackgroundTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    """Test case for the CircularApertureFlux algorithm/plugin."""

    def setUp(self):
        self.algName = "base_LocalBackground"
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0),
                                        lsst.afw.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(49.5, 49.5))

        self.bgValue = 1234.56789
        self.bgStdev = 12.3456789
        self.annulusInner = 5
        self.annulusOuter = 15
        psfSize = self.dataset.psfShape.getDeterminantRadius()
        self.numPix = np.pi*(self.annulusOuter**2 - self.annulusInner**2)*psfSize**2

    def tearDown(self):
        del self.bbox
        del self.dataset

    def checkCatalog(self, task, catalog, exposure, *args):
        exposure.maskedImage.image.array[:] += self.bgValue
        task.run(catalog, exposure, *args)
        self.assertEqual(len(catalog), 1)
        src = catalog[0]

        # Check the background measurements
        bgValue = src.get(self.algName + "_flux")
        bgStdev = src.get(self.algName + "_fluxSigma")
        self.assertFalse(src.get(self.algName + "_flag"))
        self.assertFalse(src.get(self.algName + "_flag_noGoodPixels"))
        self.assertFalse(src.get(self.algName + "_flag_noPsf"))
        self.assertFloatsAlmostEqual(bgValue, self.bgValue,
                                     atol=3.0*self.bgStdev/np.sqrt(self.numPix))  # Within 3 sigma
        self.assertFloatsAlmostEqual(bgStdev, self.bgStdev, rtol=0.1)

        # Check that the background value is useful, and it's what we'd get if the background wasn't there.
        src = catalog[0]
        psfFlux = src.get("base_PsfFlux_flux") - src.get("base_PsfFlux_area")*src.get(self.algName + "_flux")
        exposure.maskedImage.image.array[:] -= self.bgValue
        task.run(catalog, exposure, *args)
        self.assertFloatsAlmostEqual(src.get("base_PsfFlux_flux"), psfFlux, rtol=3.0e-4)

    def setConfig(self, config):
        config.plugins[self.algName].annulusInner = self.annulusInner
        config.plugins[self.algName].annulusOuter = self.annulusOuter
        config.plugins.names.add("base_PsfFlux")

    def testSingleFramePlugin(self):
        config = self.makeSingleFrameMeasurementConfig(self.algName)
        self.setConfig(config)
        task = self.makeSingleFrameMeasurementTask(config=config)
        exposure, catalog = self.dataset.realize(self.bgStdev, task.schema, randomSeed=12345)
        self.checkCatalog(task, catalog, exposure)

    def testForcedPlugin(self):
        config = self.makeForcedMeasurementConfig(self.algName)
        self.setConfig(config)
        task = self.makeForcedMeasurementTask(self.algName, config=config)
        measWcs = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs(), randomSeed=12345)
        measDataset = self.dataset.transform(measWcs)
        exposure, truthCatalog = measDataset.realize(self.bgStdev, measDataset.makeMinimalSchema(),
                                                     randomSeed=12345)
        refCat = self.dataset.catalog
        refWcs = self.dataset.exposure.getWcs()
        measCat = task.generateMeasCat(exposure, refCat, refWcs)
        task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)
        self.checkCatalog(task, measCat, exposure, refCat, refWcs)


class LocalBackgroundTransformTestCase(FluxTransformTestCase, SingleFramePluginTransformSetupHelper,
                                       lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.LocalBackgroundAlgorithm.Control
    algorithmClass = LocalBackgroundAlgorithm
    transformClass = lsst.meas.base.LocalBackgroundTransform
    flagNames = ('flag', 'flag_noGoodPixels')
    singleFramePlugins = ('base_LocalBackground',)
    forcedPlugins = ('base_LocalBackground',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
