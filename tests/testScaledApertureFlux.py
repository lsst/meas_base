#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
import math

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.base.tests

from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)

class ScaledApertureFluxTestCase(AlgorithmTestCase):
    def setUp(self):
        self.center = lsst.afw.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0),
                                        lsst.afw.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.sourceFlux = 100000.0
        self.dataset.addSource(self.sourceFlux, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def makeAlgorithm(self, ctrl=None):
        """
        Construct a ScaledApertureFluxAlgorithm with an accompanying schema and return them both.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.ScaledApertureFluxControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.ScaledApertureFluxAlgorithm(ctrl, "base_ScaledApertureFlux", schema)
        return algorithm, schema

    def testSourceFlux(self):
        """
        Check that we recover the source flux.
        """
        ctrl = lsst.meas.base.ScaledApertureFluxControl()
        algorithm, schema = self.makeAlgorithm(ctrl)
        exposure, catalog = self.dataset.realize(10.0, schema)

        # Default aperture should collect ~all source flux.
        algorithm.measure(catalog[0], exposure)
        self.assertAlmostEqual(catalog[0].get("base_ScaledApertureFlux_flux") / self.sourceFlux, 1.0, 2)
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))

        # Aperture equal to the PSF FWHM should collect ~93.7% of the flux.
        ctrl.scale = 1.0
        algorithm, schema = self.makeAlgorithm(ctrl)
        algorithm.measure(catalog[0], exposure)
        self.assertAlmostEqual(catalog[0].get("base_ScaledApertureFlux_flux") / self.sourceFlux, 0.9375, 2)
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))

    def testApertureTruncated(self):
        """
        Check that we set a flag appropriately when the aperture overflows the image.

        Note that this is a fatal failure: we do not return a useful result, but rather set the global flag.
        """
        ctrl = lsst.meas.base.ScaledApertureFluxControl()
        ctrl.scale = 100
        algorithm, schema = self.makeAlgorithm(ctrl)
        exposure, catalog = self.dataset.realize(10.0, schema)

        algorithm.measure(catalog[0], exposure)
        self.assertTrue(math.isnan(catalog[0].get("base_ScaledApertureFlux_flux")))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))

    def testSincCoeffsTruncated(self):
        """
        Check that we set a flag appropriately when the coefficient image is clipped.

        Note that we don't regard this as a fatal failure, so the global flag
        is not set and we still provide a numeric result.
        """
        ctrl = lsst.meas.base.ScaledApertureFluxControl()
        ctrl.scale = 10
        algorithm, schema = self.makeAlgorithm(ctrl)
        exposure, catalog = self.dataset.realize(10.0, schema)

        algorithm.measure(catalog[0], exposure)
        self.assertFalse(math.isnan(catalog[0].get("base_ScaledApertureFlux_flux")))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))


class ScaledApertureFluxTransformTestCase(FluxTransformTestCase, SingleFramePluginTransformSetupHelper):
    controlClass = lsst.meas.base.ScaledApertureFluxControl
    algorithmClass = lsst.meas.base.ScaledApertureFluxAlgorithm
    transformClass = lsst.meas.base.ScaledApertureFluxTransform
    flagNames = ('flag', 'flag_apertureTruncated', 'flag_sincCoeffsTruncated')
    singleFramePlugins = ('base_ScaledApertureFlux',)
    forcedPlugins = ('base_ScaledApertureFlux',)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ScaledApertureFluxTestCase)
    suites += unittest.makeSuite(ScaledApertureFluxTransformTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
