#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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

from contextlib import contextmanager
import unittest

import lsst.geom
import lsst.daf.base
import lsst.meas.base
import lsst.utils.tests

from lsst.meas.base.tests import AlgorithmTestCase


@contextmanager
def onlyLogFatal(log):
    """
    For the duration of the context, only log FATAL errors.

    This is convenient when testing algorithms under failure conditions: we
    want to be able to check that they have set appropriate flags without
    spewing alarming & confusing error messages to the console.
    """
    oldLevel = log.getThreshold()
    log.setThreshold(log.FATAL)
    try:
        yield
    finally:
        log.setThreshold(oldLevel)


class BlendednessTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(1, 4),
                                    lsst.geom.Extent2I(110, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        with self.dataset.addBlend() as family:
            family.addChild(flux=2E5, centroid=lsst.geom.Point2D(47, 33))
            family.addChild(flux=1.5E5, centroid=lsst.geom.Point2D(53, 31))

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testAbsExpectation(self):
        f = lsst.meas.base.BlendednessAlgorithm.computeAbsExpectation
        # comparison values computed with Mathematica
        self.assertFloatsAlmostEqual(f(-1.0, 1.5**2), 0.897767011, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(0.0, 1.5**2), 1.19682684, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(1.0, 1.5**2), 1.64102639, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(-1.0, 0.3**2), 0.0783651228, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(0.0, 0.3**2), 0.239365368, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(1.0, 0.3**2), 1.00046288, rtol=1E-5)

    def testAbsBias(self):
        f = lsst.meas.base.BlendednessAlgorithm.computeAbsBias
        # comparison values computed with Mathematica
        self.assertFloatsAlmostEqual(f(0.0, 1.5**2), 1.19682684, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(0.5, 1.5**2), 0.762708343, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(1.0, 1.5**2), 0.453358941, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(0.0, 0.3**2), 0.239365368, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(0.5, 0.3**2), 0.011895931, rtol=1E-5)
        self.assertFloatsAlmostEqual(f(1.0, 0.3**2), 0.0000672467314, rtol=1E-5)

    def testBlendedness(self):
        """
        Check that we measure a positive blendedness for two overlapping sources
        """
        task = self.makeSingleFrameMeasurementTask("base_Blendedness")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        self.assertGreater(catalog[1].get('base_Blendedness_abs_flux'), 0)
        self.assertGreater(catalog[2].get('base_Blendedness_abs_flux'), 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
