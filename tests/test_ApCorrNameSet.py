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

import lsst.utils.tests
import lsst.meas.base


class MinimalistTestPlugin(lsst.meas.base.SingleFramePlugin):
    """This class is used below to test registration. It is set up as the
    minimal class that still has a valid implementation. Whilst these
    methods are not needed in this test file, the class registered here
    might appear in other tests that scan the registry."""
    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def fail(self, measRecord, error=None):
        pass


class ApCorrNameTestCase(lsst.utils.tests.TestCase):

    def testDefaultNames(self):
        apCorrSet = lsst.meas.base.getApCorrNameSet()
        self.assertIn("base_PsfFlux", apCorrSet)
        self.assertIn("base_GaussianFlux", apCorrSet)

    def testAdd(self):
        nameSet0 = lsst.meas.base.getApCorrNameSet()

        lsst.meas.base.addApCorrName("test_NewName")
        nameSet1 = lsst.meas.base.getApCorrNameSet()
        self.assertIn("test_NewName", nameSet1)
        self.assertEqual(len(nameSet1 - nameSet0), 1)

        # adding a name twice is silently ignored
        lsst.meas.base.addApCorrName("test_NewName")
        nameSet2 = lsst.meas.base.getApCorrNameSet()
        self.assertIn("test_NewName", nameSet2)
        self.assertEqual(len(nameSet2 - nameSet1), 0)

    def testCopy(self):
        """Make sure getApCorrNameSet returns a copy
        """
        nameSet0 = lsst.meas.base.getApCorrNameSet()
        nameSet0.add("test_LocalName")
        nameSet1 = lsst.meas.base.getApCorrNameSet()
        self.assertNotIn("test_LocalName", nameSet1)
        self.assertEqual(len(nameSet0 - nameSet1), 1)

    def testRegisterDecorator(self):
        """Test the shouldApCorr argument of the register decorator for measurement plugins."""
        @lsst.meas.base.register("test_ApCorrPlugin", shouldApCorr=True)
        class ApCorrPlugin(MinimalistTestPlugin):
            pass

        @lsst.meas.base.register("test_NonApCorrPlugin")
        class NonApCorrPlugin(MinimalistTestPlugin):
            pass

        apCorrSet = lsst.meas.base.getApCorrNameSet()
        self.assertIn("test_ApCorrPlugin", apCorrSet)
        self.assertNotIn("test_NonApCorrPlugin", apCorrSet)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
