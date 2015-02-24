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

import numpy
import unittest
import numpy

import lsst.utils.tests
import lsst.meas.base.tests

class SdssShapeTestCase(lsst.meas.base.tests.AlgorithmTestCase):

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(240, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def assertFinite(self, value):
        self.assertTrue(numpy.isfinite(value), msg="%s is not finite" % value)

    def testGaussians(self):
        """Test that we get correct shape when measuring Gaussians with known position."""
        task = self.makeSingleFrameMeasurementTask("base_SdssShape")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(exposure, catalog)
        key = lsst.meas.base.SdssShapeResultKey(catalog.schema["base_SdssShape"])
        for record in catalog:
            result = record.get(key)
            self.assertClose(result.flux, record.get("truth_flux"), rtol=1E-2)
            self.assertFinite(result.fluxSigma)
            self.assertClose(result.x, record.get("truth_x"), rtol=1E-2)
            self.assertClose(result.y, record.get("truth_y"), rtol=1E-2)
            self.assertClose(result.xx, record.get("truth_xx"), rtol=1E-2)
            self.assertClose(result.yy, record.get("truth_yy"), rtol=1E-2)
            self.assertClose(result.xy, record.get("truth_xy"), rtol=1E-2, atol=2E-2)
            self.assertFinite(result.xxSigma)
            self.assertFinite(result.yySigma)
            self.assertFinite(result.xySigma)
            self.assertFinite(result.flux_xx_Cov)
            self.assertFinite(result.flux_yy_Cov)
            self.assertFinite(result.flux_xy_Cov)
            self.assertFinite(result.xx_yy_Cov)
            self.assertFinite(result.xx_xy_Cov)
            self.assertFinite(result.yy_xy_Cov)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SdssShapeTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
