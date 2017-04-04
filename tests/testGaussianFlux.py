#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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

from __future__ import division, absolute_import, print_function
from builtins import zip
from builtins import range
import unittest

import numpy as np

import lsst.meas.base
import lsst.utils.tests
from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)


class GaussianFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(240, 1600))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def makeAlgorithm(self, ctrl=None):
        """Construct an algorithm (finishing a schema in the process), and return both."""
        if ctrl is None:
            ctrl = lsst.meas.base.GaussianFluxControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.GaussianFluxAlgorithm(ctrl, "base_GaussianFlux", schema)
        return algorithm, schema

    def testGaussians(self):
        """Test that we get correct fluxes when measuring Gaussians with known positions and shapes."""
        task = self.makeSingleFrameMeasurementTask("base_GaussianFlux")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(catalog, exposure)
        for measRecord in catalog:
            self.assertClose(measRecord.get("base_GaussianFlux_flux"),
                             measRecord.get("truth_flux"), rtol=3E-3)

    def testMonteCarlo(self):
        """Test that we get exactly the right answer on an ideal sim with no noise, and that
        the reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        algorithm, schema = self.makeAlgorithm()
        exposure, catalog = self.dataset.realize(1E-8, schema)
        record = catalog[0]
        flux = record.get("truth_flux")
        algorithm.measure(record, exposure)
        self.assertClose(record.get("base_GaussianFlux_flux"), flux, rtol=1E-3)
        self.assertLess(record.get("base_GaussianFlux_fluxSigma"), 1E-3)
        for noise in (0.001, 0.01, 0.1):
            fluxes = []
            fluxSigmas = []
            nSamples = 1000
            for repeat in range(nSamples):
                exposure, catalog = self.dataset.realize(noise*flux, schema)
                record = catalog[1]
                algorithm.measure(record, exposure)
                fluxes.append(record.get("base_GaussianFlux_flux"))
                fluxSigmas.append(record.get("base_GaussianFlux_fluxSigma"))
            fluxMean = np.mean(fluxes)
            fluxSigmaMean = np.mean(fluxSigmas)
            fluxStandardDeviation = np.std(fluxes)
            self.assertClose(fluxSigmaMean, fluxStandardDeviation, rtol=0.10)   # rng dependent
            self.assertLess(fluxMean - flux, 2.0*fluxSigmaMean / nSamples**0.5)   # rng dependent

    def testForcedPlugin(self):
        task = self.makeForcedMeasurementTask("base_GaussianFlux")
        measWcs = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs())
        measDataset = self.dataset.transform(measWcs)
        exposure, truthCatalog = measDataset.realize(10.0, measDataset.makeMinimalSchema())
        refWcs = self.dataset.exposure.getWcs()
        refCat = self.dataset.catalog
        measCat = task.generateMeasCat(exposure, refCat, refWcs)
        task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)
        task.run(measCat, exposure, refCat, refWcs)
        for measRecord, truthRecord in zip(measCat, truthCatalog):
            # Centroid tolerances set to ~ single precision epsilon
            self.assertClose(measRecord.get("slot_Centroid_x"), truthRecord.get("truth_x"), rtol=1E-7)
            self.assertClose(measRecord.get("slot_Centroid_y"), truthRecord.get("truth_y"), rtol=1E-7)
            self.assertFalse(measRecord.get("base_GaussianFlux_flag"))
            # GaussianFlux isn't designed to do a good job in forced mode, because it doesn't account
            # for changes in the PSF (and in fact confuses them with changes in the WCS).  Hence, this
            # is really just a regression test, with the initial threshold set to just a bit more than
            # what it was found to be at one point.
            self.assertClose(measRecord.get("base_GaussianFlux_flux"), truthCatalog.get("truth_flux"),
                             rtol=0.3)
            self.assertLess(measRecord.get("base_GaussianFlux_fluxSigma"), 500.0)


class GaussianFluxTransformTestCase(FluxTransformTestCase, SingleFramePluginTransformSetupHelper,
                                    lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.GaussianFluxControl
    algorithmClass = lsst.meas.base.GaussianFluxAlgorithm
    transformClass = lsst.meas.base.GaussianFluxTransform
    singleFramePlugins = ('base_GaussianFlux',)
    forcedPlugins = ('base_GaussianFlux',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
