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

import unittest

import lsst.utils.tests
import lsst.meas.base.tests

class ClassificationTestCase(lsst.meas.base.tests.AlgorithmTestCase):

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -20),
                                        lsst.afw.geom.Extent2I(250, 150))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def testSingleFramePlugin(self):
        config = self.makeSingleFrameMeasurementConfig("base_ClassificationExtendedness",
                                                       dependencies=["base_PsfFlux"])
        # n.b. we use the truth value as ModelFlux
        config.slots.psfFlux = "base_PsfFlux"
        task = self.makeSingleFrameMeasurementTask(config=config)
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        task.run(exposure, catalog)
        self.assertLess(catalog[0].get("base_ClassificationExtendedness_value"), 0.5)
        self.assertGreater(catalog[1].get("base_ClassificationExtendedness_value"), 0.5)

    def testFlags(self):
        """Test all the failure modes of this algorithm, as well as checking that it succeeds when it should.

        Since this algorithm depends on having a ModelFlux and a PsfFlux measurement, it is a failure
        mode when either is NAN, or when ModelFluxFlag or PsfFluxFlag is True.

        When psfFluxFactor != 0, the PsfFluxErr cannot be NAN, but otherwise is ignored

        When modelFluxFactor != 0, the ModelFluxErr cannot be NAN, but otherwise is ignored
        """
        config = self.makeSingleFrameMeasurementConfig("base_ClassificationExtendedness",
                                                       dependencies=["base_PsfFlux", "base_GaussianFlux"])
        config.slots.psfFlux = "base_PsfFlux"
        config.slots.modelFlux = "base_GaussianFlux"

        def runFlagTest(psfFlux=100.0, modelFlux=200.0,
                       psfFluxSigma=1.0, modelFluxSigma=2.0,
                       psfFluxFlag=False, modelFluxFlag=False):
            task = self.makeSingleFrameMeasurementTask(config=config)
            exposure, catalog = self.dataset.realize(10.0, task.schema)
            source = catalog[0]
            source.set("base_PsfFlux_flux", psfFlux)
            source.set("base_PsfFlux_fluxSigma", psfFluxSigma)
            source.set("base_PsfFlux_flag", psfFluxFlag)
            source.set("base_GaussianFlux_flux", modelFlux)
            source.set("base_GaussianFlux_fluxSigma", modelFluxSigma)
            source.set("base_GaussianFlux_flag", modelFluxFlag)
            task.plugins["base_ClassificationExtendedness"].measure(source, exposure)
            return source.get("base_ClassificationExtendedness_flag")

        #  Test no error case - all necessary values are set
        self.assertFalse(runFlagTest())

        #  Test psfFlux flag case - failure in PsfFlux
        self.assertTrue(runFlagTest(psfFluxFlag=True))

        #  Test modelFlux flag case - failure in ModelFlux
        self.assertTrue(runFlagTest(modelFluxFlag=True))

        #  Test modelFlux NAN case
        self.assertTrue(runFlagTest(modelFlux=float("NaN"), modelFluxFlag=True))

        #  Test psfFlux NAN case
        self.assertTrue(runFlagTest(psfFlux=float("NaN"), psfFluxFlag=True))

        #  Test modelFluxErr NAN case when modelErrFactor is zero and non-zero
        config.plugins["base_ClassificationExtendedness"].modelErrFactor = 0.
        self.assertFalse(runFlagTest(modelFluxSigma=float("NaN")))
        config.plugins["base_ClassificationExtendedness"].modelErrFactor = 1.
        self.assertTrue(runFlagTest(modelFluxSigma=float("NaN")))

        #  Test psfFluxErr NAN case when psfErrFactor is zero and non-zero
        config.plugins["base_ClassificationExtendedness"].psfErrFactor = 0.
        self.assertFalse(runFlagTest(psfFluxSigma=float("NaN")))
        config.plugins["base_ClassificationExtendedness"].psfErrFactor = 1.
        self.assertTrue(runFlagTest(psfFluxSigma=float("NaN")))

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ClassificationTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
