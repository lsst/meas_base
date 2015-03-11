#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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

import math
import os
import unittest
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.utils.tests as utilsTests

from lsst.meas.base.tests import AlgorithmTestCase, TransformTestCase

class GaussianFluxTestCase(measBase.tests.AlgorithmTestCase):
    def testAlgorithm(self):
        sfmConfig = measBase.sfm.SingleFrameMeasurementConfig()
        mapper = afwTable.SchemaMapper(self.truth.getSchema())
        mapper.addMinimalSchema(self.truth.getSchema())
        outSchema = mapper.getOutputSchema()
        #  Basic test of GaussianFlux algorithm, no C++ slots
        sfmConfig.plugins = ["base_SdssCentroid", "base_GaussianFlux", "base_SdssShape"]
        sfmConfig.slots.centroid = "base_SdssCentroid"
        sfmConfig.slots.shape = "base_SdssShape"
        sfmConfig.slots.psfFlux = None
        sfmConfig.slots.instFlux = None
        sfmConfig.slots.apFlux = None
        sfmConfig.slots.modelFlux = "base_GaussianFlux"
        task = measBase.sfm.SingleFrameMeasurementTask(outSchema, config=sfmConfig)
        measCat = afwTable.SourceCatalog(outSchema)
        measCat.extend(self.truth, mapper=mapper)
        # now run the SFM task with the test plugin
        task.run(measCat, self.calexp)
        for record in measCat[:2]:
            # check all the flags
            self.assertFalse(record.get("base_GaussianFlux_flag"))
            # check the slots
            flux = record.get("base_GaussianFlux_flux")
            fluxerr = record.get("base_GaussianFlux_fluxSigma")
            truthFlux = record.get("truth_flux")
            # if a star, see if the flux measured is decent
            if record.get("truth_isStar"):
                self.assertClose(truthFlux, flux, atol=None, rtol=.02)
            if (not record.get("base_GaussianFlux_flag")):
                self.assertEqual(record.getModelFlux(), flux)
                self.assertEqual(record.getModelFluxErr(), fluxerr)
        # on the third test source (which is a blended parent), the SdssShape algorithm fails because the
        # centroid moves around too much, which gives us an opportunity to test GaussianFlux's error handling
        # note that because the shape measurement is made now even if the shape flag is set, the global
        # error flag should not be set
        for record in measCat[2:3]:
            self.assertFalse(record.get("base_GaussianFlux_flag"))
            self.assertTrue(record.get("base_GaussianFlux_flag_badShape"))
            self.assertFalse(record.get("base_GaussianFlux_flag_badCentroid"))


class GaussianFluxTransformTestCase(TransformTestCase):
    controlClass = measBase.GaussianFluxControl
    algorithmClass = measBase.GaussianFluxAlgorithm
    transformClass = measBase.GaussianFluxTransform

    def testTransform(self):
        flux, fluxSigma = 10, 1 # Arbitrary values for testing
        for flag in (True, False):
            record = self.inputCat.addNew()
            record[self.name + '_flux'] = flux
            record[self.name + '_fluxSigma'] = fluxSigma
            record.set(self.name + '_flag', flag)

        self._runTransform()

        # We are not testing the Calib itself; we assume that it is correct
        # when converting flux to magnitude, and merely check that the
        # transform has used it properly
        mag, magErr = self.calexp.getCalib().getMagnitude(flux, fluxSigma)
        for inSrc, outSrc in zip(self.inputCat, self.outputCat):
            self.assertEqual(outSrc[self.name + '_mag'], mag)
            self.assertEqual(outSrc[self.name + '_magErr'], magErr)
            self.assertEqual(outSrc.get(self.name + '_flag'), inSrc.get(self.name + '_flag'))


def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(GaussianFluxTestCase)
    suites += unittest.makeSuite(GaussianFluxTransformTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
