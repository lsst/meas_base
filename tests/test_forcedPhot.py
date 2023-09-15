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

"""Tests of the various forced photometry tasks.

These tests primarily confirm that their respective Tasks can be configured and
run without errors, but do not check anything about their algorithmic quality.
"""

import unittest

import numpy as np

import lsst.afw.image
from lsst.afw.math import ChebyshevBoundedField
from lsst.afw.table import CoordKey
from lsst.meas.base import ForcedPhotCcdTask, ForcedPhotCcdFromDataFrameTask
import lsst.meas.base.tests
import lsst.utils.tests

skyCenter = lsst.geom.SpherePoint(245.0, -45.0, lsst.geom.degrees)


class ForcedPhotometryTests:
    """Base class for tests of forced photometry tasks.

    Creates a simple test image and catalog to run forced photometry on.
    """
    def setUp(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox, crval=skyCenter)
        dataset.addSource(instFlux=1000, centroid=lsst.geom.Point2D(30, 30))
        dataset.addSource(instFlux=10000, centroid=lsst.geom.Point2D(60, 70))

        schema = dataset.makeMinimalSchema()
        # Add coordinate error fields needed by updateSourceCoords below:
        CoordKey.addErrorFields(schema)
        self.exposure, self.refCat = dataset.realize(noise=10, schema=schema)
        lsst.afw.table.updateSourceCoords(self.exposure.wcs, self.refCat)
        # Simple aperture correction map in case the task needs it.
        apCorrMap = lsst.afw.image.ApCorrMap()
        apCorrMap["base_PsfFlux_instFlux"] = ChebyshevBoundedField(bbox, np.array([[2.0]]))
        apCorrMap["base_PsfFlux_instFluxErr"] = ChebyshevBoundedField(bbox, np.array([[3.0]]))
        self.exposure.info.setApCorrMap(apCorrMap)

        # Offset WCS so that the forced coordinates don't match the truth.
        self.offsetWcs = dataset.makePerturbedWcs(self.exposure.wcs)


class ForcedPhotCcdTaskTestCase(ForcedPhotometryTests, lsst.utils.tests.TestCase):
    def testRun(self):
        config = ForcedPhotCcdTask.ConfigClass()
        task = ForcedPhotCcdTask(refSchema=self.refCat.schema, config=config)
        measCat = task.measurement.generateMeasCat(self.exposure, self.refCat, self.exposure.wcs)

        task.run(measCat, self.exposure, self.refCat, self.offsetWcs)

        # Check that something was measured.
        self.assertTrue(np.isfinite(measCat["base_TransformedCentroid_x"]).all())
        self.assertTrue(np.isfinite(measCat["base_TransformedCentroid_y"]).all())
        self.assertTrue(np.isfinite(measCat["base_PsfFlux_instFlux"]).all())
        # We use an offset WCS, so the transformed centroids should not exactly
        # match the original positions.
        self.assertFloatsNotEqual(measCat["base_TransformedCentroid_x"], self.refCat['truth_x'])
        self.assertFloatsNotEqual(measCat["base_TransformedCentroid_y"], self.refCat['truth_y'])


class ForcedPhotCcdFromDataFrameTaskTestCase(ForcedPhotometryTests, lsst.utils.tests.TestCase):
    def testRun(self):
        """Testing run() for this task ignores the dataframe->SourceCatalog
        conversion that happens in runQuantum, but that should be tested
        separately.
        """
        config = ForcedPhotCcdFromDataFrameTask.ConfigClass()
        task = ForcedPhotCcdFromDataFrameTask(refSchema=self.refCat.schema, config=config)
        measCat = task.measurement.generateMeasCat(self.exposure, self.refCat, self.exposure.wcs)

        task.run(measCat, self.exposure, self.refCat, self.offsetWcs)

        # Check that something was measured.
        self.assertTrue(np.isfinite(measCat["base_TransformedCentroidFromCoord_x"]).all())
        self.assertTrue(np.isfinite(measCat["base_TransformedCentroidFromCoord_y"]).all())
        self.assertTrue(np.isfinite(measCat["base_PsfFlux_instFlux"]).all())
        # We use an offset WCS, so the transformed centroids should not exactly
        # match the original positions.
        self.assertFloatsNotEqual(measCat["base_TransformedCentroidFromCoord_x"], self.refCat['truth_x'])
        self.assertFloatsNotEqual(measCat["base_TransformedCentroidFromCoord_y"], self.refCat['truth_y'])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
