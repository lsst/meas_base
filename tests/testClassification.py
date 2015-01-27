#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 LSST Corporation.
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
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementTask
from lsst.meas.base.base import *
from lsst.meas.base.tests import *
import unittest
import lsst.utils.tests
import numpy

numpy.random.seed(1234)


DATA_DIR = os.path.join(os.environ["MEAS_BASE_DIR"], "tests")

class SFMTestCase(lsst.meas.base.tests.AlgorithmTestCase):

    def testAlgorithm(self):

        srccat, bbox = MakeTestData.makeCatalog()
        exposure = MakeTestData.makeEmptyExposure(bbox)
        MakeTestData.fillImages(srccat, exposure)

        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        mapper = SchemaMapper(srccat.getSchema())
        mapper.addMinimalSchema(srccat.getSchema())
        outschema = mapper.getOutputSchema()

        #  Basic test of Classification algorithm, no C++ slots
        sfm_config.plugins.names = ["base_SdssCentroid", "base_PsfFlux", "base_CircularApertureFlux",
                                    "base_ClassificationExtendedness"]
        sfm_config.plugins["base_CircularApertureFlux"].radii = [7.0]
        sfm_config.slots.centroid = "base_SdssCentroid"
        sfm_config.slots.shape = None
        sfm_config.slots.psfFlux = "base_PsfFlux"
        sfm_config.slots.modelFlux = "base_CircularApertureFlux_0"
        sfm_config.slots.apFlux = None
        sfm_config.slots.instFlux = None
        task = SingleFrameMeasurementTask(outschema, config=sfm_config)
        measCat = SourceCatalog(outschema)
        measCat.extend(srccat, mapper=mapper)
        # now run the SFM task with the test plugin
        task.run(measCat, exposure)
        for i in range(len(measCat)):
            record = measCat[i]
            srcRec = srccat[i]
            # check all the flags
            # check the slots
            # if a star, see if the flux measured is decent
            probability = record.get("base_ClassificationExtendedness_value")
            if srcRec.get("truth_isStar"):
                self.assertEqual(probability, 0.0)
            else:
                self.assertEqual(probability, 1.0)

    def testFlags(self):
        """Test all the failure modes of this algorithm, as well as checking that it succeeds when it should.

        Since this algorithm depends on having a ModelFlux and a PsfFlux measurement, it is a failure
        mode when either is NAN, or when ModelFluxFlag or PsfFluxFlag is True.

        When psfFluxFactor != 0, the PsfFluxErr cannot be NAN, but otherwise is ignored

        When modelFluxFactor != 0, the ModelFluxErr cannot be NAN, but otherwise is ignored
        """
        exp = afwImage.ExposureF()
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()

        sfm_config = lsst.meas.base.sfm.SingleFrameMeasurementConfig()
        sfm_config.plugins.names = ["base_SdssCentroid", "base_PsfFlux", "base_CircularApertureFlux",
                                    "base_ClassificationExtendedness"]
        sfm_config.plugins["base_CircularApertureFlux"].radii = [7.0]
        sfm_config.slots.centroid = "base_SdssCentroid"
        sfm_config.slots.shape = None
        sfm_config.slots.psfFlux = "base_PsfFlux"
        sfm_config.slots.modelFlux = "base_CircularApertureFlux_0"
        sfm_config.slots.apFlux = None
        sfm_config.slots.instFlux = None
        task = SingleFrameMeasurementTask(schema, config=sfm_config)
        measCat = SourceCatalog(schema)

        #  Test no error case - all necessary values are set
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        source.set("base_CircularApertureFlux_0_fluxSigma", 2)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertFalse(source.get("base_ClassificationExtendedness_flag"))

        #  Test psfFlux flag case - failure in PsfFlux
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        source.set("base_CircularApertureFlux_0_fluxSigma", 2)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        source.set("base_PsfFlux_flag", True)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertTrue(source.get("base_ClassificationExtendedness_flag"))

        #  Test modelFlux flag case - falure in ModelFlux
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        source.set("base_CircularApertureFlux_0_fluxSigma", 2)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        source.set("base_CircularApertureFlux_0_flag", True)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertTrue(source.get("base_ClassificationExtendedness_flag"))

        #  Test modelFlux NAN case
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_fluxSigma", 2)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        source.set("base_CircularApertureFlux_0_flag", True)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertTrue(source.get("base_ClassificationExtendedness_flag"))

        #  Test psfFlux NAN case
        source = measCat.addNew()
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        source.set("base_CircularApertureFlux_0_fluxSigma", 2)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        source.set("base_CircularApertureFlux_0_flag", True)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertTrue(source.get("base_ClassificationExtendedness_flag"))

        #  Test modelFluxErr NAN case when modelErrFactor is zero and non-zero
        sfm_config.plugins["base_ClassificationExtendedness"].modelErrFactor = 0.
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertFalse(source.get("base_ClassificationExtendedness_flag"))

        sfm_config.plugins["base_ClassificationExtendedness"].modelErrFactor = 1.
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_PsfFlux_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertTrue(source.get("base_ClassificationExtendedness_flag"))

        #  Test psfFluxErr NAN case when psfErrFactor is zero and non-zero
        sfm_config.plugins["base_ClassificationExtendedness"].psfErrFactor = 0.
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_CircularApertureFlux_0_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertFalse(source.get("base_ClassificationExtendedness_flag"))

        sfm_config.plugins["base_ClassificationExtendedness"].psfErrFactor = 1.
        source = measCat.addNew()
        source.set("base_PsfFlux_flux", 100)
        source.set("base_CircularApertureFlux_0_fluxSigma", 1)
        source.set("base_CircularApertureFlux_0_flux", 200)
        task.plugins["base_ClassificationExtendedness"].measure(source, exp)
        self.assertTrue(source.get("base_ClassificationExtendedness_flag"))

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SFMTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
