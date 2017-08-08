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

from __future__ import absolute_import, division, print_function
import unittest

import lsst.meas.base
import lsst.utils.tests

from lsst.meas.base.tests import FluxTransformTestCase, SingleFramePluginTransformSetupHelper


class PeakLikelihoodFluxTransformTestCase(FluxTransformTestCase,
                                          SingleFramePluginTransformSetupHelper,
                                          lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.PeakLikelihoodFluxControl
    algorithmClass = lsst.meas.base.PeakLikelihoodFluxAlgorithm
    transformClass = lsst.meas.base.PeakLikelihoodFluxTransform
    singleFramePlugins = ('base_PeakLikelihoodFlux',)
    forcedPlugins = ('base_PeakLikelihoodFlux',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
