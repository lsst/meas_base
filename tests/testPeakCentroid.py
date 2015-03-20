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

import numpy
import unittest

import lsst.daf.base
import lsst.meas.base
import lsst.utils.tests

from lsst.meas.base.tests import (CentroidTransformTestCase, SingleFramePluginTransformSetupHelper,
                                  ForcedPluginTransformSetupHelper)

class SingleFramePeakCentroidTestCase(CentroidTransformTestCase, SingleFramePluginTransformSetupHelper):
    class SingleFramePeakCentroidPluginFactory(object):
        """
        Helper class to sub in an empty PropertyList as the final argument to
        lsst.meas.base.SingleFramePeakCentroidPlugin.
        """
        def __call__(self, control, name, inputSchema):
            return lsst.meas.base.SingleFramePeakCentroidPlugin(control, name, inputSchema,
                                                                lsst.daf.base.PropertyList())
    controlClass = lsst.meas.base.SingleFramePeakCentroidConfig
    algorithmClass = SingleFramePeakCentroidPluginFactory()
    transformClass = lsst.meas.base.SimpleCentroidTransform
    flagNames = ()
    singleFramePlugins = ("base_PeakCentroid",)


class ForcedPeakCentroidTestCase(CentroidTransformTestCase, ForcedPluginTransformSetupHelper):
    controlClass = lsst.meas.base.ForcedPeakCentroidConfig
    algorithmClass = lsst.meas.base.ForcedPeakCentroidPlugin
    transformClass = lsst.meas.base.SimpleCentroidTransform
    flagNames = ()
    forcedPlugins = ("base_PeakCentroid",)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SingleFramePeakCentroidTestCase)
    suites += unittest.makeSuite(ForcedPeakCentroidTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
