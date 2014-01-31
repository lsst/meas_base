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

import os
import unittest
import numpy

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.base

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data")

class PsfFluxTestCase(lsst.utils.tests.TestCase):
    """Test case for the PsfFlux algorithm/plugins
    """

    def setUp(self):
        self.truth = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-01.fits"))
        self.calexp = lsst.afw.image.ExposureF.readFits(os.path.join(DATA_DIR, "calexp-01.fits"))

    def tearDown(self):
        del self.truth
        del self.calexp

    def testDirect(self):
        """Verify that we get reasonable answers when we call apply() directly on an isolated source,
        and that these results are identical regardless of which overload we use.
        """
        record = self.truth[0]
        centroidKey = self.truth.schema.find("truth.centroid").key
        result1 = lsst.meas.base.PsfFluxAlgorithm.apply(
            self.calexp, record.getFootprint(),
            record.get(centroidKey)
            )
        result2 = lsst.meas.base.PsfFluxAlgorithm.apply(
            self.calexp,
            lsst.meas.base.PsfFluxAlgorithm.Input(record.getFootprint(), record.get(centroidKey))
            )
        self.assertEqual(result1.flux, result2.flux)
        self.assertEqual(result1.fluxSigma, result2.fluxSigma)
        self.assertFalse(result1.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_PSF))
        self.assertFalse(result1.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_GOOD_PIXELS))
        self.assertFalse(result1.getFlag(lsst.meas.base.PsfFluxAlgorithm.EDGE))
        self.assertFalse(result2.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_PSF))
        self.assertFalse(result2.getFlag(lsst.meas.base.PsfFluxAlgorithm.NO_GOOD_PIXELS))
        self.assertFalse(result2.getFlag(lsst.meas.base.PsfFluxAlgorithm.EDGE))



def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
