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
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.utils.tests


class CoordErrorConventionTestCase(lsst.utils.tests.TestCase):
    """Verify that ``coord_raErr`` / ``coord_decErr`` / ``coord_ra_dec_Cov``
    carry the tangent-plane uncertainty propagated through a local gnomonic
    Jacobian.
    """

    def setUp(self):
        # Reference position at high declination so cos(Dec) = 0.5
        self.referenceRa = 0.0 * lsst.geom.degrees
        self.referenceDec = 60.0 * lsst.geom.degrees
        self.cosDec = np.cos(self.referenceDec.asRadians())

        self.pixelScale = 0.2 * lsst.geom.arcseconds
        self.pixelScaleRadians = self.pixelScale.asRadians()

        crpix = lsst.geom.Point2D(0.0, 0.0)
        crval = lsst.geom.SpherePoint(self.referenceRa, self.referenceDec)
        cdMatrix = afwGeom.makeCdMatrix(scale=self.pixelScale,
                                        orientation=0.0 * lsst.geom.radians,
                                        flipX=False)
        self.wcs = afwGeom.makeSkyWcs(crpix, crval, cdMatrix)

        # Minimal SourceCatalog with a centroid (slot + uncertainty fields)
        # and the global coord_raErr / coord_decErr / coord_ra_dec_Cov fields
        # that SourceRecord.updateCoord writes into.
        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("test_x", type="D", units="pixel",
                        doc="centroid x")
        schema.addField("test_y", type="D", units="pixel",
                        doc="centroid y")
        schema.addField("test_xErr", type="F", units="pixel",
                        doc="1-sigma uncertainty on x")
        schema.addField("test_yErr", type="F", units="pixel",
                        doc="1-sigma uncertainty on y")
        schema.addField("test_x_y_Cov", type="F", units="pixel^2",
                        doc="x/y covariance")
        afwTable.CoordKey.addErrorFields(schema)
        schema.getAliasMap().set("slot_Centroid", "test")

        self.catalog = afwTable.SourceCatalog(schema)

    def _addSource(self, sigmaPix, xyCov=0.0):
        """Add one source at the WCS reference pixel with covariance
        ``[[sigmaPix**2, xyCov], [xyCov, sigmaPix**2]]``.
        """
        source = self.catalog.addNew()
        source["test_x"] = 0.0
        source["test_y"] = 0.0
        source["test_xErr"] = sigmaPix
        source["test_yErr"] = sigmaPix
        source["test_x_y_Cov"] = xyCov
        return source

    def test_skyPositionMatchesWcsReference(self):
        """Check the simple case where the recorded sky position is the WCS
        reference.
        """
        source = self._addSource(sigmaPix=1.0)
        source.updateCoord(self.wcs)
        self.assertAlmostEqual(source["coord_ra"], self.referenceRa.asRadians(),
                               places=10)
        self.assertAlmostEqual(source["coord_dec"], self.referenceDec.asRadians(),
                               places=10)

    def test_coordErrIsTangentPlaneSigma(self):
        """Check that the coordinate errors equal the pixel-scale sigma.
        """
        sigmaPix = 1.0
        source = self._addSource(sigmaPix)
        source.updateCoord(self.wcs)
        expectedTangentErr = sigmaPix * self.pixelScaleRadians
        self.assertFloatsAlmostEqual(source["coord_decErr"], expectedTangentErr, rtol=1e-3)

        # Note that there is no cos(Dec) factor.
        self.assertFloatsAlmostEqual(source["coord_raErr"], expectedTangentErr, rtol=1e-3)

    def test_raDecCovIsTangentPlaneCov(self):
        """Check the calculations of the covariance matrix.

        Note that ``coord_ra_dec_Cov`` is Cov(RA*cos(Dec), Dec).
        """
        sigmaPix = 1.0

        # Diagonal pixel covariance -> diagonal sky covariance.
        source = self._addSource(sigmaPix, xyCov=0.0)
        source.updateCoord(self.wcs)
        self.assertFloatsAlmostEqual(source["coord_ra_dec_Cov"], 0.0,
                                     atol=1e-20)

        # Non-diagonal pixel covariance.
        rhoPix = 0.5
        source = self._addSource(sigmaPix, xyCov=rhoPix * sigmaPix**2)
        source.updateCoord(self.wcs)
        expectedCov = self.pixelScaleRadians**2 * (rhoPix * sigmaPix**2)
        if not self.wcs.isFlipped:
            expectedCov *= -1
        self.assertFloatsAlmostEqual(source["coord_ra_dec_Cov"], expectedCov, rtol=1e-3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
