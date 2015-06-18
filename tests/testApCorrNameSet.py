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

import lsst.utils.tests
import lsst.meas.base

class ApCorrNameTestCase(unittest.TestCase):

    def testDefaultNames(self):
        apCorrSet = lsst.meas.base.getApCorrNameSet()
        self.assertTrue("base_PsfFlux" in apCorrSet)
        self.assertTrue("base_GaussianFlux" in apCorrSet)

    def testAdd(self):
        nameSet0 = lsst.meas.base.getApCorrNameSet()

        lsst.meas.base.addApCorrName("test_NewName")
        nameSet1 = lsst.meas.base.getApCorrNameSet()
        self.assertTrue("test_NewName" in nameSet1)
        self.assertEqual(len(nameSet1 - nameSet0), 1)

        # adding a name twice is silently ignored
        lsst.meas.base.addApCorrName("test_NewName")
        nameSet2 = lsst.meas.base.getApCorrNameSet()
        self.assertTrue("test_NewName" in nameSet2)
        self.assertEqual(len(nameSet2 - nameSet1), 0)

    def testCopy(self):
        """Make sure getApCorrNameSet returns a copy
        """
        nameSet0 = lsst.meas.base.getApCorrNameSet()
        nameSet0.add("test_LocalName")
        nameSet1 = lsst.meas.base.getApCorrNameSet()
        self.assertTrue("test_LocalName" not in nameSet1)
        self.assertEqual(len(nameSet0 - nameSet1), 1)

    def testRegisterDecorator(self):
        """Test the canApCorr argument of the register decorator for measurement plugins
        """
        @lsst.meas.base.register("test_ApCorrPlugin", canApCorr=True)
        class ApCorrPlugin(lsst.meas.base.SingleFramePlugin):
            pass

        @lsst.meas.base.register("test_NonApCorrPlugin")
        class NonApCorrPlugin(lsst.meas.base.SingleFramePlugin):
            pass

        apCorrSet = lsst.meas.base.getApCorrNameSet()
        self.assertTrue("test_ApCorrPlugin" in apCorrSet)
        self.assertFalse("test_NonApCorrPlugin" in apCorrSet)



def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ApCorrNameTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
