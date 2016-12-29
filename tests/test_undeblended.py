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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests for measuring sources on undeblended images
"""

from __future__ import absolute_import, division, print_function

import sys
import unittest

import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.base as measBase
import lsst.utils.tests

class UndeblendedTestCase(lsst.utils.tests.TestCase):
    def testUndeblendedMeasurement(self):
        """Check undeblended measurement and aperture correction"""
        width, height = 100, 100  # Dimensions of image
        x0, y0 = 1234, 5678  # Offset of image
        radius = 3.0  # Aperture radius
        xCenter, yCenter = width//2, height//2  # Position of first source; integer values, for convenience
        xOffset, yOffset = 1, 1  # Offset from first source to second source
        flux1, flux2 = 1000, 1  # Flux of sources
        apCorrValue = 3.21  # Aperture correction value to apply

        image = afwImage.MaskedImageF(afwGeom.ExtentI(width, height))
        image.setXY0(x0, y0)
        image.getVariance().set(1.0)

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("centroid_x", type=float)
        schema.addField("centroid_y", type=float)
        schema.addField("centroid_flag", type='Flag')
        schema.getAliasMap().set("slot_Centroid", "centroid")

        sfmConfig = measBase.SingleFrameMeasurementConfig()
        algName = "base_CircularApertureFlux"

        for subConfig in (sfmConfig.plugins, sfmConfig.undeblended):
            subConfig.names = [algName]
            subConfig[algName].radii = [radius]
            subConfig[algName].maxSincRadius = 0  # Disable sinc photometry because we're undersampled
        slots = sfmConfig.slots
        slots.centroid = "centroid"
        slots.shape = None
        slots.apFlux = None
        slots.modelFlux = None
        slots.psfFlux = None
        slots.instFlux = None
        slots.calibFlux = None

        fieldName = lsst.meas.base.CircularApertureFluxAlgorithm.makeFieldPrefix(algName, radius)
        measBase.addApCorrName(fieldName)

        apCorrConfig = measBase.ApplyApCorrConfig()
        apCorrConfig.proxies = {"undeblended_" + fieldName: fieldName}

        sfm = measBase.SingleFrameMeasurementTask(config=sfmConfig, schema=schema)
        apCorr = measBase.ApplyApCorrTask(config=apCorrConfig, schema=schema)

        cat = afwTable.SourceCatalog(schema)
        parent = cat.addNew()
        parent.set("centroid_x", x0 + xCenter)
        parent.set("centroid_y", y0 + yCenter)
        parent.setFootprint(afwDetection.Footprint(afwGeom.Point2I(x0 + xCenter, y0 + yCenter), radius))

        # First child is bright, dominating the blend
        child1 = cat.addNew()
        child1.set("centroid_x", parent.get("centroid_x"))
        child1.set("centroid_y", parent.get("centroid_y"))
        child1.setParent(parent.getId())
        image.set(xCenter, yCenter, (flux1, 0, 0))
        foot1 = afwDetection.Footprint(afwGeom.Point2I(x0 + xCenter, y0 + yCenter), 0.1)
        child1.setFootprint(afwDetection.HeavyFootprintF(foot1, image))

        # Second child is fainter, but we want to be able to measure it!
        child2 = cat.addNew()
        child2.set("centroid_x", parent.get("centroid_x") + xOffset)
        child2.set("centroid_y", parent.get("centroid_y") + yOffset)
        child2.setParent(parent.getId())
        image.set(xCenter + xOffset, yCenter + yOffset, (flux2, 0, 0))
        foot2 = afwDetection.Footprint(afwGeom.Point2I(x0 + xCenter + xOffset, y0 + yCenter + yOffset), 0.1)
        child2.setFootprint(afwDetection.HeavyFootprintF(foot2, image))

        spans = []
        for ss in foot1.getSpans():
            spans.append(ss)
        for ss in foot2.getSpans():
            spans.append(ss)
        bbox = afwGeom.Box2I()
        bbox.include(foot1.getBBox())
        bbox.include(foot2.getBBox())
        parent.setFootprint(afwDetection.Footprint(spans, bbox))

        exposure = afwImage.makeExposure(image)

        sfm.run(cat, exposure)

        def checkSource(source, baseName, expectedFlux):
            """Check that we get the expected results"""
            self.assertEqual(source.get(baseName + "_flux"), expectedFlux)
            self.assertGreater(source.get(baseName + "_fluxSigma"), 0)
            self.assertFalse(source.get(baseName + "_flag"))

        # Deblended
        checkSource(child1, fieldName, flux1)
        checkSource(child2, fieldName, flux2)

        # Undeblended
        checkSource(child1, "undeblended_" + fieldName, flux1 + flux2)
        checkSource(child2, "undeblended_" + fieldName, flux1 + flux2)

        # Apply aperture correction
        apCorrMap = afwImage.ApCorrMap()
        apCorrMap[fieldName + "_flux"] = afwMath.ChebyshevBoundedField(
            image.getBBox(),
            apCorrValue*np.ones((1, 1), dtype=float)
        )
        apCorrMap[fieldName + "_fluxSigma"] = afwMath.ChebyshevBoundedField(
            image.getBBox(),
            apCorrValue*np.zeros((1, 1), dtype=float)
        )

        apCorr.run(cat, apCorrMap)

        # Deblended
        checkSource(child1, fieldName, flux1*apCorrValue)
        checkSource(child2, fieldName, flux2*apCorrValue)

        # Undeblended
        checkSource(child1, "undeblended_" + fieldName, (flux1 + flux2)*apCorrValue)
        checkSource(child2, "undeblended_" + fieldName, (flux1 + flux2)*apCorrValue)

        self.assertIn(fieldName + "_apCorr", schema)
        self.assertIn(fieldName + "_apCorrSigma", schema)
        self.assertIn("undeblended_" + fieldName + "_apCorr", schema)
        self.assertIn("undeblended_" + fieldName + "_apCorrSigma", schema)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    setup_module(sys.modules[__name__])
    unittest.main()
