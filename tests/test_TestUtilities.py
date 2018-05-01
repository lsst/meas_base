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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#


import numpy as np
import unittest
import lsst.afw.geom as afwGeom
import lsst.utils.tests as utilsTests

# Rename the class on import so it does not confuse the test scanner
from lsst.meas.base.tests import TestDataset as DatasetTester


class TestDatasetTestCase(utilsTests.TestCase):

    def setUp(self):
        # Construct an arbitrary WCS for testing.
        crval = afwGeom.SpherePoint(45.0, 45.0, afwGeom.degrees)
        scale = 0.2*afwGeom.arcseconds
        crpix = afwGeom.PointD(100, 100)
        self.wcs = afwGeom.makeSkyWcs(crpix=crpix, crval=crval,
                                      cdMatrix=afwGeom.makeCdMatrix(scale=scale))

    def tearDown(self):
        del self.wcs

    def test_perturb(self):
        """Test that perturbing a WCS gives us back something different."""
        # We should always get something different from our starting point.
        self.assertNotEqual(self.wcs, DatasetTester.makePerturbedWcs(self.wcs))

        # If we use the same random seed, the results should be reproducible.
        self.assertEqual(DatasetTester.makePerturbedWcs(self.wcs, randomSeed=0),
                         DatasetTester.makePerturbedWcs(self.wcs, randomSeed=0))

        # If we specify different seeds, we should always get something
        # different.
        self.assertNotEqual(DatasetTester.makePerturbedWcs(self.wcs, randomSeed=1),
                            DatasetTester.makePerturbedWcs(self.wcs, randomSeed=2))

    def test_randomState(self):
        """Test that we do not alter global state when perturbing the WCS."""
        # This checks that global state doesn't change while the test is
        # executing. This isn't perfectly robust: the test will fail if
        # another thread manipulates the RNG state while this is executing.
        # However, since pytest doesn't using multi-threading, it should be
        # safe in practice.
        init_state = np.random.get_state()
        DatasetTester.makePerturbedWcs(self.wcs)
        for init, final in zip(init_state, np.random.get_state()):
            self.assertTrue(np.array_equal(init, final))


class TestMemory(utilsTests.MemoryTestCase):
    pass


def setup_module(module):
    utilsTests.init()


if __name__ == "__main__":
    utilsTests.init()
    unittest.main()
