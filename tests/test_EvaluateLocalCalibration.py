# This file is part of ap_association.
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

import numpy as np
import unittest

from lsst.afw.cameraGeom.testUtils import DetectorWrapper
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.geom as geom
import lsst.utils.tests

from lsst.meas.base import EvaluateLocalCalibrationTask


class TestLocalPhotoCalibration(lsst.meas.base.tests.AlgorithmTestCase,
                                lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testNoFlags(self):
        task = self.makeSingleFrameMeasurementTask("base_LocalPhotoCalib")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]

        calib = exposure.getPhotoCalib().getLocalCalibration(self.center)
        calibErr = exposure.getPhotoCalib().getCalibrationErr()
        self.assertEqual(record.get("base_LocalPhotoCalib"), calib)
        self.assertEqual(record.get("base_LocalPhotoCalibErr"), calibErr)


class TestLocalWcs(lsst.meas.base.tests.AlgorithmTestCase,
                   lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def testNoFlags(self):
        task = self.makeSingleFrameMeasurementTask("base_LocalWcs")
        exposure, catalog = self.dataset.realize(10.0,
                                                 task.schema,
                                                 randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]

        localCDMatrix = exposure.getWcs().getCdMatrix(self.center)
        self.assertEqual(record.get("base_LocalWcs_CDMatrix_1_1"),
                         localCDMatrix[0, 0])
        self.assertEqual(record.get("base_LocalWcs_CDMatrix_2_1"),
                         localCDMatrix[1, 0])
        self.assertEqual(record.get("base_LocalWcs_CDMatrix_1_2"),
                         localCDMatrix[0, 1])
        self.assertEqual(record.get("base_LocalWcs_CDMatrix_2_2"),
                         localCDMatrix[1, 1])


class TestEvaluateLocalCalibrationTask(unittest.TestCase):

    def setUp(self):
        """Create a sqlite3 database with default tables and schemas.
        """
        # metadata taken from CFHT data
        # v695856-e0/v695856-e0-c000-a00.sci_img.fits
        afwImageUtils.resetFilters()
        afwImageUtils.defineFilter('g', lambdaEff=487, alias="g.MP9401")

        self.metadata = dafBase.PropertySet()

        self.metadata.set("SIMPLE", "T")
        self.metadata.set("BITPIX", -32)
        self.metadata.set("NAXIS", 2)
        self.metadata.set("NAXIS1", 1024)
        self.metadata.set("NAXIS2", 1153)
        self.metadata.set("RADECSYS", 'FK5')
        self.metadata.set("EQUINOX", 2000.)

        self.metadata.setDouble("CRVAL1", 215.604025685476)
        self.metadata.setDouble("CRVAL2", 53.1595451514076)
        self.metadata.setDouble("CRPIX1", 1109.99981456774)
        self.metadata.setDouble("CRPIX2", 560.018167811613)
        self.metadata.set("CTYPE1", 'RA---SIN')
        self.metadata.set("CTYPE2", 'DEC--SIN')

        self.metadata.setDouble("CD1_1", 5.10808596133527E-05)
        self.metadata.setDouble("CD1_2", 1.85579539217196E-07)
        self.metadata.setDouble("CD2_2", -5.10281493481982E-05)
        self.metadata.setDouble("CD2_1", -8.27440751733828E-07)

        self.wcs = afwGeom.makeSkyWcs(self.metadata)
        self.exposure = afwImage.makeExposure(
            afwImage.makeMaskedImageFromArrays(np.ones((1024, 1153))),
            self.wcs)
        detector = DetectorWrapper(id=23, bbox=self.exposure.getBBox()).detector
        visit = afwImage.VisitInfo(
            exposureId=1234,
            exposureTime=200.,
            date=dafBase.DateTime("2014-05-13T17:00:00.000000000",
                                  dafBase.DateTime.Timescale.TAI))
        self.exposure.setDetector(detector)
        self.exposure.getInfo().setVisitInfo(visit)
        self.exposure.setFilter(afwImage.Filter('g'))
        self.flux0 = 10000
        self.flux0_err = 100

        self.calibration = 1e-3
        self.calibrationErr = 1e-4
        self.linearXCalibration = lsst.afw.math.ChebyshevBoundedField(
            self.exposure.getBBox(),
            np.array([[self.calibration,
                       self.calibration]]))
        photoCalib = lsst.afw.image.PhotoCalib(self.linearXCalibration,
                                               self.calibrationErr)
        self.exposure.setPhotoCalib(photoCalib)

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema.addField("base_Centroid_x", type="D")
        self.schema.addField("base_Centroid_y", type="D")

    def test_addColumns(self):
        """Test that schema columns are added properly.
        """
        EvaluateLocalCalibrationTask(self.schema)

        # Find photoCalib and wcs columns in schema.
        try:
            self.schema.find("base_localPhotoCalib")
            self.schema.find("base_localPhotoCalibErr")

            self.schema.find("base_CDMatrix_1_1")
            self.schema.find("base_CDMatrix_1_2")
            self.schema.find("base_CDMatrix_2_1")
            self.schema.find("base_CDMatrix_2_2")
        except KeyError:
            self.fail("Calibration schema columns not properly added.")

    def test_run(self):
        evalLocCalib = EvaluateLocalCalibrationTask(self.schema)
        sourceCat = self._create_sources(5)

        evalLocCalib.run(sourceCat, self.exposure)

        wcs = self.exposure.getWcs()
        photoCalib = self.exposure.getPhotoCalib()

        for idx, srcRec in enumerate(sourceCat):
            localCalib = photoCalib.getLocalCalibration(
                geom.Point2D(100 * idx, 100 * idx))
            self.assertAlmostEqual(srcRec.get("base_localPhotoCalib"),
                                   localCalib)
            self.assertAlmostEqual(srcRec.get("base_localPhotoCalibErr"),
                                   self.calibrationErr)

            cdMatrix = wcs.getCdMatrix(geom.Point2D(100 * idx, 100 * idx))
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_1_1"),
                                   cdMatrix[0, 0])
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_1_2"),
                                   cdMatrix[0, 1])
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_2_1"),
                                   cdMatrix[1, 0])
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_2_2"),
                                   cdMatrix[1, 1])

    def test_wcsStorage(self):
        """Test that the correct wcs values are stored in the source catalog.
        """
        evalLocCalib = EvaluateLocalCalibrationTask(self.schema)
        sourceCat = self._create_sources(5)
        wcs = self.exposure.getWcs()

        for idx, srcRec in enumerate(sourceCat):
            evalLocCalib._storeWcs(srcRec, wcs)
            cdMatrix = wcs.getCdMatrix(geom.Point2D(100 * idx, 100 * idx))
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_1_1"),
                                   cdMatrix[0, 0])
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_1_2"),
                                   cdMatrix[0, 1])
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_2_1"),
                                   cdMatrix[1, 0])
            self.assertAlmostEqual(srcRec.get("base_CDMatrix_2_2"),
                                   cdMatrix[1, 1])

    def _create_sources(self, n_sources):
        """Create sources with centroids and ra/dec.

        Parameters
        ----------
        n_sources : `int`
            Number of sources to create.
        """
        sourceCat = afwTable.SourceCatalog(self.schema)
        for idx in range(5):
            rec = sourceCat.addNew()
            rec.set("id", idx)
            rec.set("base_Centroid_x", 100 * idx)
            rec.set("base_Centroid_y", 100 * idx)

        sourceCat.defineCentroid("base_Centroid")
        return sourceCat

    def test_photoCalibStorage(self):
        """Test that the correct PhotoCalib values are stored in the source
        catalog.
        """
        evalLocCalib = EvaluateLocalCalibrationTask(self.schema)
        sourceCat = self._create_sources(5)
        photoCalib = self.exposure.getPhotoCalib()

        for idx, srcRec in enumerate(sourceCat):
            evalLocCalib._storePhotoCalib(srcRec, photoCalib)
            localCalib = photoCalib.getLocalCalibration(
                geom.Point2D(100 * idx, 100 * idx))
            self.assertAlmostEqual(srcRec.get("base_localPhotoCalib"),
                                   localCalib)
            self.assertAlmostEqual(srcRec.get("base_localPhotoCalibErr"),
                                   self.calibrationErr)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
