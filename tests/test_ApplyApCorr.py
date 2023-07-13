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

import lsst.utils.tests
import lsst.geom
import lsst.meas.base.tests
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.base.applyApCorr as applyApCorr
from lsst.afw.math import ChebyshevBoundedField
from lsst.meas.base.apCorrRegistry import addApCorrName


def initializeSourceCatalog(schema=None, name=None, instFlux=None, sigma=None, centroid=None):
    instFluxName = name + "_instFlux"
    instFluxErrName = name + "_instFluxErr"
    instFluxKey = schema.find(instFluxName).key
    centroidKey = afwTable.Point2DKey(schema["slot_Centroid"])
    sourceCat = afwTable.SourceCatalog(schema)
    source = sourceCat.addNew()
    source.set(instFluxKey, instFlux)
    source.set(instFluxErrName, sigma)
    source.set(centroidKey, centroid)
    return sourceCat


class ApplyApCorrTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        names = ["test2", "test"]
        for name in names:
            addApCorrName(name)
            schema.addField(name + "_instFlux", type=np.float64)
            schema.addField(name + "_instFluxErr", type=np.float64)
            schema.addField(name + "_flag", type=np.float64)
            schema.addField(name + "_Centroid_x", type=np.float64)
            schema.addField(name + "_Centroid_y", type=np.float64)
            schema.getAliasMap().set('slot_Centroid', name + '_Centroid')
        self.ap_corr_task = applyApCorr.ApplyApCorrTask(schema=schema)
        self.name = name   # just use 'test' prefix for most of the tests
        self.schema = schema

    def tearDown(self):
        del self.schema
        del self.ap_corr_task

    def testAddFields(self):
        # Check that the required fields have been added to the schema
        self.assertIn(self.name + "_apCorr", self.schema.getNames())
        self.assertIn(self.name + "_apCorrErr", self.schema.getNames())
        self.assertIn(self.name + "_flag_apCorr", self.schema.getNames())
        self.assertLess(self.schema.find("test_apCorr").key.getOffset(),
                        self.schema.find("test2_apCorr").key.getOffset())

    def testSuccessUnflagged(self):
        # Check that the aperture correction flag is set to False if aperture
        # correction was successfully run
        flagName = self.name + "_flag_apCorr"
        flagKey = self.schema.find(flagName).key
        source_test_instFlux = 5.1
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, instFlux=source_test_instFlux,
                                            sigma=0, centroid=source_test_centroid)
        instFluxName = self.name + "_instFlux"
        instFluxErrName = self.name + "_instFluxErr"

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10))
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[instFluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[instFluxErrName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)
        self.assertFalse(sourceCat[flagKey])

    def testFailureFlagged(self):
        # Check that aperture correction flag is set to True if aperture
        # correction is invalid (negative)
        flagName = self.name + "_flag_apCorr"
        flagKey = self.schema.find(flagName).key
        source_test_instFlux = 5.2
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, instFlux=source_test_instFlux,
                                            sigma=0, centroid=source_test_centroid)
        instFluxName = self.name + "_instFlux"
        instFluxErrName = self.name + "_instFluxErr"

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10))
        coefficients = -(np.ones((1, 1), dtype=np.float64))
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[instFluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[instFluxErrName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)
        self.assertTrue(sourceCat[flagKey])

    def testCatFluxUnchanged(self):
        # Pick arbitrary but unique values for the test case
        source_test_instFlux = 5.3
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, instFlux=source_test_instFlux,
                                            sigma=0, centroid=source_test_centroid)
        instFluxName = self.name + "_instFlux"
        instFluxErrName = self.name + "_instFluxErr"
        instFluxKey = self.schema.find(instFluxName).key

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10))
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[instFluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[instFluxErrName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)

        self.assertEqual(sourceCat[instFluxKey], source_test_instFlux)

    def testCatFluxHalf(self):
        # Pick arbitrary but unique values for the test case
        source_test_instFlux = 5.4
        source_test_centroid = lsst.geom.Point2D(5, 7.1)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, instFlux=source_test_instFlux,
                                            sigma=0, centroid=source_test_centroid)
        instFluxName = self.name + "_instFlux"
        instFluxErrName = self.name + "_instFluxErr"
        instFluxKey = self.schema.find(instFluxName).key

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10))
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients /= 2.
        coefficients_sigma = np.zeros((1, 1), dtype=np.float64)
        apCorrMap[instFluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[instFluxErrName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)

        self.assertAlmostEqual(sourceCat[instFluxKey], source_test_instFlux / 2)

    def testCatFluxErr(self):
        """Test catalog flux errors.

        Notes
        -----
        This test will break if ``UseNaiveFluxErr = False``~

        The alternate method significantly overestimates noise, causing this
        test to fail.  It is likely that this test will need to be modified if
        the noise calculation is updated.
        """
        # Pick arbitrary but unique values for the test case
        source_test_instFlux = 5.5
        source_test_sigma = 0.23
        source_test_centroid = lsst.geom.Point2D(5, 7.3)
        sourceCat = initializeSourceCatalog(schema=self.schema, name=self.name, instFlux=source_test_instFlux,
                                            sigma=source_test_sigma, centroid=source_test_centroid)

        instFluxName = self.name + "_instFlux"
        instFluxErrName = self.name + "_instFluxErr"
        instFluxErrKey = self.schema.find(instFluxErrName).key

        apCorrMap = afwImage.ApCorrMap()
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.ExtentI(10, 10))
        coefficients = np.ones((1, 1), dtype=np.float64)
        coefficients_sigma = np.ones((1, 1), dtype=np.float64)
        apCorrMap[instFluxName] = ChebyshevBoundedField(bbox, coefficients)
        apCorrMap[instFluxErrName] = ChebyshevBoundedField(bbox, coefficients_sigma)
        self.ap_corr_task.run(sourceCat, apCorrMap)

        self.assertAlmostEqual(sourceCat[instFluxErrKey], source_test_sigma)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
