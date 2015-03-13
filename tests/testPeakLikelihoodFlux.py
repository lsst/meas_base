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

import lsst.meas.base
import lsst.utils.tests

from lsst.meas.base.tests import TransformTestCase

class PeakLikelihoodFluxTransformTestCase(TransformTestCase):
    controlClass = lsst.meas.base.PeakLikelihoodFluxControl
    algorithmClass = lsst.meas.base.PeakLikelihoodFluxAlgorithm
    transformClass = lsst.meas.base.PeakLikelihoodFluxTransform
    singleFramePlugins = ('base_PeakLikelihoodFlux',)
    forcedPlugins = ('base_PeakLikelihoodFlux',)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PeakLikelihoodFluxTransformTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
