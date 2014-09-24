#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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
import numpy

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.base.tests

class ApertureFluxTestCase(lsst.utils.tests.TestCase):
    """Test case for the ApertureFlux algorithm/plugins
    """

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(20, -100), lsst.afw.geom.Point2I(100, -20))
        self.exposure = lsst.afw.image.ExposureF(self.bbox)
        self.exposure.getMaskedImage().getImage().set(1.0)
        self.exposure.getMaskedImage().getVariance().set(0.25)
        self.ctrl = lsst.meas.base.ApertureFluxAlgorithm.Control()

    def tearDown(self):
        del self.bbox
        del self.exposure

    def computeNaiveArea(self, position, radius):
        """Return the area of a circular aperture with the given position and radius, according to
        the 'naive' definition of the aperture - just test whether the center of each pixel is within
        the circle.
        """
        x, y = numpy.meshgrid(numpy.arange(self.bbox.getBeginX(), self.bbox.getEndX()),
                              numpy.arange(self.bbox.getBeginY(), self.bbox.getEndY()))
        return ((x - position.getX())**2 + (y - position.getY())**2 <= radius**2).sum()

    def testNaive(self):
        positions = [lsst.afw.geom.Point2D(60.0, -60.0),
                     lsst.afw.geom.Point2D(60.5, -60.0),
                     lsst.afw.geom.Point2D(60.0, -60.5),
                     lsst.afw.geom.Point2D(60.5, -60.5)]
        radii = [12.0, 17.0]
        maskedImage = self.exposure.getMaskedImage()
        image = maskedImage.getImage()
        for position in positions:
            for radius in radii:
                ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(radius, radius, 0.0),
                                                         position)
                area = self.computeNaiveArea(position, radius)
                self.assertClose(
                    lsst.meas.base.ApertureFluxAlgorithm.computeNaiveFlux(image, ellipse, self.ctrl),
                    area
                )
                self.assertClose(
                    lsst.meas.base.ApertureFluxAlgorithm.computeFlux(image, ellipse, self.ctrl),
                    area
                )

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ApertureFluxTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
