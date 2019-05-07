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
import math

import lsst.geom
import lsst.afw.image
import lsst.afw.table
from lsst.meas.base import SincCoeffsD
from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)
import lsst.utils.tests


@unittest.skipIf(SincCoeffsD.DISABLED_AT_COMPILE_TIME, "Sinc photometry is disabled.")
class ScaledApertureFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                    lsst.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.sourceFlux = 100000.0
        self.dataset.addSource(self.sourceFlux, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def makeAlgorithm(self, ctrl=None):
        """Construct an and return both it and its schema.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.ScaledApertureFluxControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.ScaledApertureFluxAlgorithm(ctrl, "base_ScaledApertureFlux", schema)
        return algorithm, schema

    def testSourceFlux(self):
        """Check that we recover the source instFlux.
        """
        ctrl = lsst.meas.base.ScaledApertureFluxControl()
        algorithm, schema = self.makeAlgorithm(ctrl)
        exposure, catalog = self.dataset.realize(10.0, schema, randomSeed=0)

        # Default aperture should collect ~all source instFlux.
        algorithm.measure(catalog[0], exposure)
        self.assertAlmostEqual(catalog[0].get("base_ScaledApertureFlux_instFlux") / self.sourceFlux, 1.0, 2)
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))

        # Aperture equal to the PSF FWHM should collect ~93.7% of the
        # instFlux.
        ctrl.scale = 1.0
        algorithm, schema = self.makeAlgorithm(ctrl)
        algorithm.measure(catalog[0], exposure)
        self.assertAlmostEqual(catalog[0].get("base_ScaledApertureFlux_instFlux") /
                               self.sourceFlux, 0.9375, 2)
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))

    def testApertureTruncated(self):
        """Check that a flag is set when the aperture overflows the image.

        Notes
        -----
        This is a fatal failure: we do not return a useful result, but rather
        set the global flag.
        """
        ctrl = lsst.meas.base.ScaledApertureFluxControl()
        ctrl.scale = 100
        algorithm, schema = self.makeAlgorithm(ctrl)
        exposure, catalog = self.dataset.realize(10.0, schema, randomSeed=1)

        algorithm.measure(catalog[0], exposure)
        self.assertTrue(math.isnan(catalog[0].get("base_ScaledApertureFlux_instFlux")))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))

    def testSincCoeffsTruncated(self):
        """Check that a flag is set when the coefficient image is clipped.

        Notes
        -----
        This is not regarded as a fatal failure, so the global flag is not set
        and we still provide a numeric result.
        """
        ctrl = lsst.meas.base.ScaledApertureFluxControl()
        ctrl.scale = 10
        algorithm, schema = self.makeAlgorithm(ctrl)
        exposure, catalog = self.dataset.realize(10.0, schema, randomSeed=2)

        algorithm.measure(catalog[0], exposure)
        self.assertFalse(math.isnan(catalog[0].get("base_ScaledApertureFlux_instFlux")))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag"))
        self.assertFalse(catalog[0].get("base_ScaledApertureFlux_flag_apertureTruncated"))
        self.assertTrue(catalog[0].get("base_ScaledApertureFlux_flag_sincCoeffsTruncated"))


class ScaledApertureFluxTransformTestCase(FluxTransformTestCase,
                                          SingleFramePluginTransformSetupHelper,
                                          lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.ScaledApertureFluxControl
    algorithmClass = lsst.meas.base.ScaledApertureFluxAlgorithm
    transformClass = lsst.meas.base.ScaledApertureFluxTransform
    flagNames = ('flag', 'flag_apertureTruncated', 'flag_sincCoeffsTruncated')
    singleFramePlugins = ('base_ScaledApertureFlux',)
    forcedPlugins = ('base_ScaledApertureFlux',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
