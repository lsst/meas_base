#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#

import unittest

import numpy

import lsst.utils.tests
import lsst.meas.base.tests
from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  CentroidTransformTestCase, SingleFramePluginTransformSetupHelper)

class SdssShapeTestCase(AlgorithmTestCase):

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(240, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))
        self.config = self.makeSingleFrameMeasurementConfig("base_SdssShape")

    def tearDown(self):
        del self.bbox
        del self.dataset
        del self.config

    def assertFinite(self, value):
        self.assertTrue(numpy.isfinite(value), msg="%s is not finite" % value)

    def _runMeasurementTask(self):
        task = self.makeSingleFrameMeasurementTask("base_SdssShape", config=self.config)
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(exposure, catalog)
        return exposure, catalog

    def _checkShape(self, result, record):
        self.assertClose(result.x, record.get("truth_x"), rtol=1E-2)
        self.assertClose(result.y, record.get("truth_y"), rtol=1E-2)
        self.assertClose(result.xx, record.get("truth_xx"), rtol=1E-2)
        self.assertClose(result.yy, record.get("truth_yy"), rtol=1E-2)
        self.assertClose(result.xy, record.get("truth_xy"), rtol=1E-1, atol=2E-1)
        self.assertFinite(result.xxSigma)
        self.assertFinite(result.yySigma)
        self.assertFinite(result.xySigma)
        self.assertFinite(result.flux_xx_Cov)
        self.assertFinite(result.flux_yy_Cov)
        self.assertFinite(result.flux_xy_Cov)
        self.assertFinite(result.xx_yy_Cov)
        self.assertFinite(result.xx_xy_Cov)
        self.assertFinite(result.yy_xy_Cov)
        self.assertFalse(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.FAILURE))
        self.assertFalse(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.UNWEIGHTED_BAD))
        self.assertFalse(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.UNWEIGHTED))
        self.assertFalse(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.SHIFT))
        self.assertFalse(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.MAXITER))

    def _checkPsfShape(self, result, psfResult, psfTruth):
        self.assertClose(psfResult.getIxx(), psfTruth.getIxx(), rtol=1E-4)
        self.assertClose(psfResult.getIyy(), psfTruth.getIyy(), rtol=1E-4)
        self.assertClose(psfResult.getIxy(), psfTruth.getIxy(), rtol=1E-4)
        self.assertFalse(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.PSF_SHAPE_BAD))

    def testMeasureGoodPsf(self):
        """Test that we measure shapes and record the PSF shape correctly."""
        exposure, catalog = self._runMeasurementTask()
        key = lsst.meas.base.SdssShapeResultKey(catalog.schema["base_SdssShape"])
        psfTruth = exposure.getPsf().computeShape()
        for record in catalog:
            result = record.get(key)
            self._checkShape(result, record)
            psfResult = key.getPsfShape(record)
            self._checkPsfShape(result, psfResult, psfTruth)

    def testMeasureWithoutPsf(self):
        """Test that we measure shapes correctly and do not record the PSF shape when not needed."""
        self.config.plugins["base_SdssShape"].doMeasurePsf = False
        _, catalog = self._runMeasurementTask()
        key = lsst.meas.base.SdssShapeResultKey(catalog.schema["base_SdssShape"])
        for record in catalog:
            result = record.get(key)
            self._checkShape(result, record)
        self.assertFalse("base_SdssShape_psf_xx" in catalog.schema)
        self.assertFalse("base_SdssShape_psf_yy" in catalog.schema)
        self.assertFalse("base_SdssShape_psf_xy" in catalog.schema)
        self.assertFalse("base_SdssShape_flag_psf" in catalog.schema)

    def testMeasureBadPsf(self):
        """Test that we measure shapes correctly and set a flag with the PSF is unavailable."""
        self.config.plugins["base_SdssShape"].doMeasurePsf = True
        task = self.makeSingleFrameMeasurementTask("base_SdssShape", config=self.config)
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        exposure.setPsf(None) # Set Psf to None to test no psf case
        task.run(exposure, catalog)
        key = lsst.meas.base.SdssShapeResultKey(catalog.schema["base_SdssShape"])
        for record in catalog:
            result = record.get(key)
            self._checkShape(result, record)
            self.assertTrue(result.getFlag(lsst.meas.base.SdssShapeAlgorithm.PSF_SHAPE_BAD))


class SdssShapeTransformTestCase(FluxTransformTestCase, CentroidTransformTestCase,
                                 SingleFramePluginTransformSetupHelper):
    name = "sdssShape"
    controlClass = lsst.meas.base.SdssShapeControl
    algorithmClass = lsst.meas.base.SdssShapeAlgorithm
    transformClass = lsst.meas.base.SdssShapeTransform
    flagNames = ("flag", "flag_unweighted", "flag_unweightedBad", "flag_shift", "flag_maxIter")
    singleFramePlugins = ('base_SdssShape',)
    forcedPlugins = ('base_SdssShape',)
    testPsf = True

    def _setFieldsInRecords(self, records, name):
        FluxTransformTestCase._setFieldsInRecords(self, records, name)
        CentroidTransformTestCase._setFieldsInRecords(self, records, name)
        for record in records:
            for field in ('xx', 'yy', 'xy', 'xxSigma', 'yySigma', 'xySigma', 'psf_xx', 'psf_yy', 'psf_xy'):
                if record.schema.join(name, field) in record.schema:
                    record[record.schema.join(name, field)] = numpy.random.random()

    def _compareFieldsInRecords(self, inSrc, outSrc, name):
        FluxTransformTestCase._compareFieldsInRecords(self, inSrc, outSrc, name)
        CentroidTransformTestCase._compareFieldsInRecords(self, inSrc, outSrc, name)

        inShape = lsst.meas.base.ShapeResultKey(inSrc.schema[name]).get(inSrc)
        outShape = lsst.meas.base.ShapeResultKey(outSrc.schema[name]).get(outSrc)

        centroid = lsst.meas.base.CentroidResultKey(inSrc.schema[name]).get(inSrc).getCentroid()
        xform = self.calexp.getWcs().linearizePixelToSky(centroid, lsst.afw.geom.radians)

        trInShape = inShape.getShape().transform(xform.getLinear())
        self.assertEqual(trInShape.getIxx(), outShape.getShape().getIxx())
        self.assertEqual(trInShape.getIyy(), outShape.getShape().getIyy())
        self.assertEqual(trInShape.getIxy(), outShape.getShape().getIxy())

        m = lsst.meas.base.makeShapeTransformMatrix(xform.getLinear())
        numpy.testing.assert_array_almost_equal(
            numpy.dot(numpy.dot(m, inShape.getShapeErr()), m.transpose()), outShape.getShapeErr()
        )

        if self.testPsf:
            inPsfShape = lsst.meas.base.ShapeResultKey(inSrc.schema[inSrc.schema.join(name, "psf")]).get(inSrc)
            outPsfShape = lsst.meas.base.ShapeResultKey(outSrc.schema[outSrc.schema.join(name, "psf")]).get(outSrc)
            trInPsfShape = inPsfShape.getShape().transform(xform.getLinear())

            self.assertEqual(trInPsfShape.getIxx(), outPsfShape.getShape().getIxx())
            self.assertEqual(trInPsfShape.getIyy(), outPsfShape.getShape().getIyy())
            self.assertEqual(trInPsfShape.getIxy(), outPsfShape.getShape().getIxy())
        else:
            self.assertNotIn(inSrc.schema.join(name, "psf", "xx"), inSrc.schema)
            self.assertNotIn(inSrc.schema.join(name, "flag", "psf"), inSrc.schema)

class PsfSdssShapeTransformTestCase(SdssShapeTransformTestCase):
    testPsf = False
    def makeSdssShapeControl(self):
        ctrl = lsst.meas.base.SdssShapeControl()
        ctrl.doMeasurePsf = False
        return ctrl
    controlClass = makeSdssShapeControl

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SdssShapeTestCase)
    suites += unittest.makeSuite(SdssShapeTransformTestCase)
    suites += unittest.makeSuite(PsfSdssShapeTransformTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
