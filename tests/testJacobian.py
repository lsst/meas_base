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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#


import unittest

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.utils.tests
import lsst.meas.base.tests

from lsst.meas.base.tests import AlgorithmTestCase
from lsst.meas.base.sfm import SingleFrameMeasurementConfig

class JacobianTestCase(AlgorithmTestCase):
    def setUp(self):
        # Pick arbitrary numbers to create a detector object, and a synthetic
        # dataset. The particular numbers have no special meaning to the test
        # and be anything as long as they are self consistent (i.e. the
        # fake source is inside the bounding box)
        self.center = lsst.afw.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(140,160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

        md = dafBase.PropertyList()

        for k, v in (
                ("EQUINOX", 2000.0),
                ("CRPIX1" , 5353.0),
                ("CRPIX2" , -35.0),
                ("CD1_1"  , 0.0),
                ("CD1_2"  , -5.611E-05),
                ("CD2_1"  , -5.611E-05),
                ("CD2_2"  , -0.0),
                ("CRVAL1" , 4.5789875),
                ("CRVAL2" , 16.30004444),
                ("CUNIT1" , 'deg'),
                ("CUNIT2" , 'deg'),
                ("CTYPE1" , 'RA---TAN'),
                ("CTYPE2" , 'DEC--TAN'),
                ):
            md.set(k, v)

        self.wcs = afwImage.makeWcs(md)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset
        del self.wcs

    def testJacobianPlugin(self):
        config = SingleFrameMeasurementConfig()
        config.plugins.names |= ['base_Jacobian']
        # Pixel scale chosen to approximately match the scale defined in the above WCS
        config.plugins['base_Jacobian'].pixelScale = 0.2
        task = self.makeSingleFrameMeasurementTask(config=config)
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        exposure.setWcs(self.wcs)
        task.run(exposure, catalog)
        record = catalog[0]
        self.assertFalse(record.get("base_Jacobian_flag"))
        self.assertAlmostEqual(record.get("base_Jacobian_value"), 1.0200183929088285, 4)

def suite():
    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(JacobianTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    "Run the tests"
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
