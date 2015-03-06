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

import unittest

import lsst.meas.base as measBase
import lsst.utils.tests as utilsTests

from lsst.meas.base.tests import TransformTestCase

class PeakLikelihoodFluxTransformTestCase(TransformTestCase):
    controlClass = measBase.PeakLikelihoodFluxControl
    algorithmClass = measBase.PeakLikelihoodFluxAlgorithm
    transformClass = measBase.PeakLikelihoodFluxTransform

    def testTransform(self):
        flux, fluxSigma = 10, 1 # Arbitrary values for testing

        for flag in (True, False):
            record = self.inputCat.addNew()
            record[self.name + '_flux'] = flux
            record[self.name + '_fluxSigma'] = fluxSigma
            record.set(self.name + '_flag', flag)

        self._runTransform()

        # We are not testing the Calib itself; we assume that it is correct
        # when converting flux to magnitude, and merely check that the
        # transform has used it properly
        mag, magErr = self.calexp.getCalib().getMagnitude(flux, fluxSigma)
        for inSrc, outSrc in zip(self.inputCat, self.outputCat):
            self.assertEqual(outSrc[self.name + '_mag'], mag)
            self.assertEqual(outSrc[self.name + '_magErr'], magErr)
            self.assertEqual(outSrc.get(self.name + '_flag'), inSrc.get(self.name + '_flag'))

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PeakLikelihoodFluxTransformTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
