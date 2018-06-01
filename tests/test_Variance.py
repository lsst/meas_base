#
# LSST Data Management System
# Copyright 2008-2017 LSST/AURA
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

"""
Tests for Variance measurement algorithm
"""
import unittest

import numpy as np

import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.meas.base as measBase
import lsst.utils.tests

try:
    display
except NameError:
    display = False


class VarianceTest(lsst.utils.tests.TestCase):

    def setUp(self):
        size = 128  # size of image (pixels)
        center = lsst.geom.Point2D(size//2, size//2)  # object center
        width = 2  # PSF width
        flux = 10.0  # Flux of object
        variance = 1.0  # Mean variance value
        varianceStd = 0.1  # Standard deviation of the variance value

        # Set a seed for predictable randomness
        np.random.seed(300)

        # Create a random image to be used as variance plane
        variancePlane = np.random.normal(variance, varianceStd, size*size).reshape(size, size)

        # Initial setup of an image
        exp = afwImage.ExposureF(size, size)
        image = exp.getMaskedImage().getImage()
        mask = exp.getMaskedImage().getMask()
        var = exp.getMaskedImage().getVariance()
        image.set(0.0)
        mask.set(0)
        var.getArray()[:, :] = variancePlane

        # Put down a PSF
        psfSize = int(6*width + 1)  # Size of PSF image; must be odd
        psf = afwDetection.GaussianPsf(psfSize, psfSize, width)
        exp.setPsf(psf)
        psfImage = psf.computeImage(center).convertF()
        psfImage *= flux
        image.Factory(image, psfImage.getBBox(afwImage.PARENT)).__iadd__(psfImage)
        var.Factory(var, psfImage.getBBox(afwImage.PARENT)).__iadd__(psfImage)

        # Put in some bad pixels to ensure they're ignored
        for i in range(-5, 6):
            bad = size//2 + i*width
            var.getArray()[bad, :] = float("nan")
            mask.getArray()[bad, :] = mask.getPlaneBitMask("BAD")
            var.getArray()[:, bad] = float("nan")
            mask.getArray()[:, bad] = mask.getPlaneBitMask("BAD")

        # Put in some unmasked bad pixels outside the expected aperture, to ensure the aperture is working
        var.getArray()[0, 0] = float("nan")
        var.getArray()[0, -1] = float("nan")
        var.getArray()[-1, 0] = float("nan")
        var.getArray()[-1, -1] = float("nan")

        if display:
            import lsst.afw.display as afwDisplay
            afwDisplay.getDisplay(1).mtv(image)
            afwDisplay.getDisplay(2).mtv(mask)
            afwDisplay.getDisplay(3).mtv(var)

        config = measBase.SingleFrameMeasurementConfig()
        config.plugins.names = ["base_NaiveCentroid", "base_SdssShape", "base_Variance"]
        config.slots.centroid = "base_NaiveCentroid"
        config.slots.psfFlux = None
        config.slots.apFlux = None
        config.slots.modelFlux = None
        config.slots.instFlux = None
        config.slots.calibFlux = None
        config.slots.shape = "base_SdssShape"
        config.slots.psfShape = None
        config.plugins["base_Variance"].mask = ["BAD", "SAT"]

        config.validate()
        schema = afwTable.SourceTable.makeMinimalSchema()

        task = measBase.SingleFrameMeasurementTask(schema, config=config)
        catalog = afwTable.SourceCatalog(schema)

        spans = afwGeom.SpanSet.fromShape(int(width))
        spans = spans.shiftedBy(int(center.getX()), int(center.getY()))
        foot = afwDetection.Footprint(spans)
        peak = foot.getPeaks().addNew()
        peak.setIx(int(center.getX()))
        peak.setIy(int(center.getY()))
        peak.setFx(center.getX())
        peak.setFy(center.getY())
        peak.setPeakValue(flux)

        source = catalog.addNew()
        source.setFootprint(foot)

        self.variance = variance
        self.varianceStd = varianceStd
        self.mask = mask
        self.catalog = catalog
        self.exp = exp
        self.task = task
        self.source = source

    def tearDown(self):
        del self.mask
        del self.catalog
        del self.exp
        del self.task
        del self.source

    def testVariance(self):
        self.task.run(self.catalog, self.exp)

        self.assertLess(np.abs(self.source.get("base_Variance_value") - self.variance), self.varianceStd)

        # flag_emptyFootprint should not have been set since the footprint has non-masked pixels at this
        # point.
        self.assertFalse(self.source.get("base_Variance_flag_emptyFootprint"))

    def testEmptyFootprint(self):
        # Set the pixel mask for all pixels to 'BAD' and remeasure.
        self.mask.getArray()[:, :] = self.mask.getPlaneBitMask("BAD")
        self.task.run(self.catalog, self.exp)

        # The computed variance should be nan and flag_emptyFootprint should have been set since the footprint
        # has all masked pixels at this point.
        self.assertTrue(np.isnan(self.source.get("base_Variance_value")))
        self.assertTrue(self.source.get("base_Variance_flag_emptyFootprint"))


class BadCentroidTest(lsst.utils.tests.TestCase):

    def testBadCentroid(self):
        """The flag from the centroid slot should propagate to the badCentroid flag on the variance plugin."""
        schema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFramePeakCentroidPlugin(measBase.SingleFramePeakCentroidConfig(),
                                               "centroid", schema, None)
        schema.getAliasMap().set("slot_Centroid", "centroid")
        variance = measBase.SingleFrameVariancePlugin(measBase.VarianceConfig(),
                                                      "variance", schema, None)
        catalog = afwTable.SourceCatalog(schema)

        # The centroid is not flagged as bad, but there's no way the algorithm can run without
        # valid data in the SourceRecord and Exposure: this should throw a logic error.
        record = catalog.addNew()
        record.set("centroid_flag", False)
        with self.assertRaises(measBase.MeasurementError) as measErr:
            variance.measure(record, None)
        variance.fail(record, measErr.exception)
        self.assertTrue(record.get("variance_flag"))
        self.assertFalse(record.get("variance_flag_badCentroid"))

        # The centroid is flagged as bad, so we should get a MeasurementError
        # indicating an expected failure.
        record = catalog.addNew()
        record.set("centroid_flag", True)
        with self.assertRaises(measBase.MeasurementError) as measErr:
            variance.measure(record, None)
        variance.fail(record, measErr.exception)
        self.assertTrue(record.get("variance_flag"))
        self.assertTrue(record.get("variance_flag_badCentroid"))

##############################################################################################################


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
