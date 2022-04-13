# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import numpy as np

import lsst.geom
import lsst.afw.geom
import lsst.meas.base
import lsst.utils.tests


class ShapeTransformMatrixTestCase(lsst.utils.tests.TestCase):

    def testIdentity(self):
        # A no-op coordinate transform translates to a no-op shape transform
        a = lsst.geom.AffineTransform()
        np.testing.assert_array_equal(a.getMatrix(), np.identity(3))
        m = lsst.meas.base.makeShapeTransformMatrix(a.getLinear())
        np.testing.assert_array_equal(m, np.identity(3))

    def testVsTransform(self):
        # Transforming an ellipse by multiplying by the matrix should be
        # equivalent to calling its transform() method.
        lt = lsst.geom.LinearTransform.makeRotation(lsst.geom.Angle(np.random.random()))
        e = lsst.afw.geom.Quadrupole(np.random.random(), np.random.random(),
                                     np.random.random())
        np.testing.assert_array_almost_equal(np.dot(lsst.meas.base.makeShapeTransformMatrix(lt),
                                                    e.getParameterVector()),
                                             e.transform(lt).getParameterVector())

    def testVales(self):
        # Test that the analytically-derived correct values are computed
        lt = lsst.geom.LinearTransform(np.random.random((2, 2)))
        m = lsst.meas.base.makeShapeTransformMatrix(lt)

        self.assertEqual(m[0, 0], lt[0, 0]*lt[0, 0])
        self.assertEqual(m[0, 1], lt[0, 1]*lt[0, 1])
        self.assertEqual(m[0, 2], 2*lt[0, 0]*lt[0, 1])

        self.assertEqual(m[1, 0], lt[1, 0]*lt[1, 0])
        self.assertEqual(m[1, 1], lt[1, 1]*lt[1, 1])
        self.assertEqual(m[1, 2], 2*lt[1, 0]*lt[1, 1])

        self.assertEqual(m[2, 0], lt[0, 0]*lt[1, 0])
        self.assertEqual(m[2, 1], lt[0, 1]*lt[1, 1])
        self.assertAlmostEqual(m[2, 2], lt[0, 0]*lt[1, 1] + lt[0, 1]*lt[1, 0], places=14)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
