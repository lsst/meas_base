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
import lsst.afw.image
import lsst.afw.table
import lsst.utils.tests

from lsst.meas.base.tests import (AlgorithmTestCase, FluxTransformTestCase,
                                  SingleFramePluginTransformSetupHelper)


def compute_chi2(exposure, centroid, instFlux, maskPlane=None):
    """Return the chi2 for the exposure PSF at the given position.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to calculate the PSF chi2 on.
    centroid : `lsst::geom::Point2D`
        Center of the source on the exposure to calculate the PSF at.
    instFlux : `float`
        Total flux of this source, as computed by the fitting algorithm.
    maskPlane : `lsst.afw.image.Image`, optional
        Pixels to mask out of the calculation.

    Returns
    -------
    chi2, nPixels: `float`, `float`
        The computed chi2 and number of pixels included in the calculation.
    """
    psfImage = exposure.getPsf().computeImage(centroid)
    scaledPsf = psfImage.array*instFlux
    # Get a sub-image the same size as the returned PSF image.
    sub = exposure.Factory(exposure, psfImage.getBBox(), lsst.afw.image.LOCAL)
    if maskPlane is not None:
        unmasked = np.logical_not(sub.mask.array & sub.mask.getPlaneBitMask(maskPlane))
        chi2 = np.sum((sub.image.array[unmasked] - scaledPsf[unmasked])**2
                      / sub.variance.array[unmasked])
        nPixels = unmasked.sum()
    else:
        chi2 = np.sum((sub.image.array - scaledPsf)**2 / sub.variance.array)
        nPixels = np.prod(sub.mask.array.shape)
    return chi2, nPixels


class PsfFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                    lsst.geom.Extent2I(100, 100))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, self.center)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.dataset

    def makeAlgorithm(self, ctrl=None):
        """Construct an algorithm and return both it and its schema.
        """
        if ctrl is None:
            ctrl = lsst.meas.base.PsfFluxControl()
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        algorithm = lsst.meas.base.PsfFluxAlgorithm(ctrl, "base_PsfFlux", schema)
        return algorithm, schema

    def testMasking(self):
        algorithm, schema = self.makeAlgorithm()
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(10.0, schema, randomSeed=0)
        record = catalog[0]
        badPoint = lsst.geom.Point2I(self.center) + lsst.geom.Extent2I(3, 4)
        imageArray = exposure.image.array
        maskArray = exposure.mask.array
        badMask = exposure.mask.getPlaneBitMask("BAD")
        imageArray[badPoint.getY() - exposure.getY0(), badPoint.getX() - exposure.getX0()] = np.inf
        maskArray[badPoint.getY() - exposure.getY0(), badPoint.getX() - exposure.getX0()] |= badMask
        # Should get an infinite value exception, because we didn't mask that
        # one pixel
        with self.assertRaises(lsst.meas.base.PixelValueError):
            algorithm.measure(record, exposure)
        # If we do mask it, we should get a reasonable result
        ctrl = lsst.meas.base.PsfFluxControl()
        ctrl.badMaskPlanes = ["BAD"]
        algorithm, schema = self.makeAlgorithm(ctrl)
        algorithm.measure(record, exposure)
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_instFlux"),
                                     record.get("truth_instFlux"),
                                     atol=3*record.get("base_PsfFlux_instFluxErr"))
        chi2, nPixels = compute_chi2(exposure, record.getCentroid(),
                                     record.get("base_PsfFlux_instFlux"),
                                     "BAD")
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_chi2"), chi2, rtol=1e-7)
        # If we mask the whole image, we should get a failure flag.
        maskArray[:, :] |= badMask
        algorithm.measure(record, exposure)
        self.assertEqual(record.get("base_PsfFlux_flag_noGoodPixels"), 1)
        self.assertEqual(record.get("base_PsfFlux_npixels"), nPixels)

    def testSubImage(self):
        """Test measurement on sub-images.

        Specifically, checks that we don't get confused by images with nonzero
        ``xy0``, and that the ``EDGE`` flag is set when it should be.
        """

        algorithm, schema = self.makeAlgorithm()
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(10.0, schema, randomSeed=1)
        record = catalog[0]
        psfImage = exposure.getPsf().computeImage(record.getCentroid())
        bbox = psfImage.getBBox()
        bbox.grow(-1)
        subExposure = exposure.Factory(exposure, bbox, lsst.afw.image.LOCAL)

        algorithm.measure(record, subExposure)
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_instFlux"), record.get("truth_instFlux"),
                                     atol=3*record.get("base_PsfFlux_instFluxErr"))
        self.assertTrue(record.get("base_PsfFlux_flag_edge"))

        # Calculating chi2 requires trimming the PSF image by one pixel per side
        # to match the subExposure created above, so we can't use compute_chi2()
        # directly here.
        scaledPsf = psfImage.array[1:-1, 1:-1]*record.get("base_PsfFlux_instFlux")
        chi2 = np.sum((subExposure.image.array - scaledPsf)**2 / subExposure.variance.array)
        nPixels = np.prod(subExposure.image.array.shape)
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_chi2"), chi2, rtol=1e-7)
        self.assertEqual(record.get("base_PsfFlux_npixels"), nPixels)

    def testNoPsf(self):
        """Test that we raise `FatalAlgorithmError` when there's no PSF.
        """
        algorithm, schema = self.makeAlgorithm()
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(10.0, schema, randomSeed=2)
        exposure.setPsf(None)
        with self.assertRaises(lsst.meas.base.FatalAlgorithmError):
            algorithm.measure(catalog[0], exposure)

    def testMonteCarlo(self):
        """Test an ideal simulation, with no noise.

        Demonstrate that:

        - We get exactly the right answer, and
        - The reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        algorithm, schema = self.makeAlgorithm()
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(0.0, schema, randomSeed=3)
        record = catalog[0]
        instFlux = record.get("truth_instFlux")
        algorithm.measure(record, exposure)
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_instFlux"), instFlux, rtol=1E-3)
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_instFluxErr"), 0.0, rtol=1E-3)
        # no noise, so infinite chi2 on this one
        self.assertEqual(record.get("base_PsfFlux_chi2"), np.inf)
        for noise in (0.001, 0.01, 0.1):
            nSamples = 1000
            instFluxes = np.zeros(nSamples, dtype=np.float64)
            instFluxErrs = np.zeros(nSamples, dtype=np.float64)
            expectedChi2s = np.zeros(nSamples, dtype=np.float64)
            expectedNPixels = -100
            measuredChi2s = np.zeros(nSamples, dtype=np.float64)
            measuredNPixels = np.zeros(nSamples, dtype=np.float64)
            for i in range(nSamples):
                # By using ``i`` to seed the RNG, we get results which
                # fall within the tolerances defined below. If we allow this
                # test to be truly random, passing becomes RNG-dependent.
                exposure, catalog = self.dataset.realize(noise*instFlux, schema, randomSeed=i)
                record = catalog[0]
                algorithm.measure(record, exposure)
                instFluxes[i] = record.get("base_PsfFlux_instFlux")
                instFluxErrs[i] = record.get("base_PsfFlux_instFluxErr")
                measuredChi2s[i] = record.get("base_PsfFlux_chi2")
                measuredNPixels[i] = record.get("base_PsfFlux_npixels")
                chi2, nPixels = compute_chi2(exposure, record.getCentroid(), instFluxes[i])
                expectedChi2s[i] = chi2
                expectedNPixels = nPixels  # should be the same each time
            instFluxMean = np.mean(instFluxes)
            instFluxErrMean = np.mean(instFluxErrs)
            instFluxStandardDeviation = np.std(instFluxes)
            self.assertFloatsAlmostEqual(instFluxErrMean, instFluxStandardDeviation, rtol=0.10)
            self.assertLess(abs(instFluxMean - instFlux), 2.0*instFluxErrMean / nSamples**0.5)
            self.assertFloatsAlmostEqual(measuredChi2s, expectedChi2s, rtol=1e-7)
            # Should have the exact same number of pixels used every time.
            self.assertTrue(np.all(measuredNPixels == expectedNPixels))

    def testSingleFramePlugin(self):
        task = self.makeSingleFrameMeasurementTask("base_PsfFlux")
        # Results are RNG dependent; we choose a seed that is known to pass.
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=4)
        task.run(catalog, exposure)
        record = catalog[0]
        self.assertFalse(record.get("base_PsfFlux_flag"))
        self.assertFalse(record.get("base_PsfFlux_flag_noGoodPixels"))
        self.assertFalse(record.get("base_PsfFlux_flag_edge"))
        self.assertFloatsAlmostEqual(record.get("base_PsfFlux_instFlux"), record.get("truth_instFlux"),
                                     atol=3*record.get("base_PsfFlux_instFluxErr"))

    def testForcedPlugin(self):
        task = self.makeForcedMeasurementTask("base_PsfFlux")
        # Results of this test are RNG dependent: we choose seeds that are
        # known to pass.
        measWcs = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs(), randomSeed=5)
        measDataset = self.dataset.transform(measWcs)
        exposure, truthCatalog = measDataset.realize(10.0, measDataset.makeMinimalSchema(), randomSeed=5)
        refCat = self.dataset.catalog
        refWcs = self.dataset.exposure.getWcs()
        measCat = task.generateMeasCat(exposure, refCat, refWcs)
        task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)
        task.run(measCat, exposure, refCat, refWcs)
        measRecord = measCat[0]
        truthRecord = truthCatalog[0]
        # Centroid tolerances set to ~ single precision epsilon
        self.assertFloatsAlmostEqual(measRecord.get("slot_Centroid_x"),
                                     truthRecord.get("truth_x"), rtol=1E-7)
        self.assertFloatsAlmostEqual(measRecord.get("slot_Centroid_y"),
                                     truthRecord.get("truth_y"), rtol=1E-7)
        self.assertFalse(measRecord.get("base_PsfFlux_flag"))
        self.assertFalse(measRecord.get("base_PsfFlux_flag_noGoodPixels"))
        self.assertFalse(measRecord.get("base_PsfFlux_flag_edge"))
        self.assertFloatsAlmostEqual(measRecord.get("base_PsfFlux_instFlux"),
                                     truthCatalog.get("truth_instFlux"), rtol=1E-3)
        self.assertLess(measRecord.get("base_PsfFlux_instFluxErr"), 500.0)


class PsfFluxTransformTestCase(FluxTransformTestCase, SingleFramePluginTransformSetupHelper,
                               lsst.utils.tests.TestCase):
    controlClass = lsst.meas.base.PsfFluxControl
    algorithmClass = lsst.meas.base.PsfFluxAlgorithm
    transformClass = lsst.meas.base.PsfFluxTransform
    flagNames = ('flag', 'flag_noGoodPixels', 'flag_edge')
    singleFramePlugins = ('base_PsfFlux',)
    forcedPlugins = ('base_PsfFlux',)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
