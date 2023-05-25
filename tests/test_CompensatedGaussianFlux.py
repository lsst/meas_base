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
import lsst.afw.geom
import lsst.meas.base
import lsst.utils.tests
from lsst.meas.base.tests import AlgorithmTestCase
from lsst.meas.base._measBaseLib import _compensatedGaussianFiltInnerProduct


class CompensatedGaussianFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    def setUp(self):
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(1000, 1000))

        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.psf_size = 2.0
        self.dataset.psfShape = lsst.afw.geom.Quadrupole(self.psf_size**2., self.psf_size**2., 0.0)

        # We want a set of point sources at various flux levels.
        self.dataset.addSource(10000.0, lsst.geom.Point2D(50.1, 49.8))
        self.dataset.addSource(20000.0, lsst.geom.Point2D(100.5, 100.4))
        self.dataset.addSource(30000.0, lsst.geom.Point2D(150.4, 149.6))
        self.dataset.addSource(40000.0, lsst.geom.Point2D(200.2, 200.3))
        self.dataset.addSource(50000.0, lsst.geom.Point2D(250.3, 250.1))
        self.dataset.addSource(60000.0, lsst.geom.Point2D(300.4, 300.2))
        self.dataset.addSource(70000.0, lsst.geom.Point2D(350.5, 350.6))
        self.dataset.addSource(80000.0, lsst.geom.Point2D(400.6, 400.0))
        self.dataset.addSource(90000.0, lsst.geom.Point2D(450.0, 450.0))
        self.dataset.addSource(100000.0, lsst.geom.Point2D(500.7, 500.8))

        # Small test for Monte Carlo
        self.bbox_single = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                           lsst.geom.Extent2I(101, 101))
        self.dataset_single = lsst.meas.base.tests.TestDataset(self.bbox_single)
        self.dataset_single.psfShape = self.dataset.psfShape

        self.dataset_single.addSource(100000.0, lsst.geom.Point2D(50.0, 50.0))

        self.all_widths = (2, 3, 5)
        self.larger_widths = (3, 5)
        self.all_ts = (1.5, 2.0)

    def tearDown(self):
        del self.bbox
        del self.dataset
        del self.bbox_single
        del self.dataset_single

    def makeAlgorithm(self, config=None):
        """Construct an algorithm and return both it and its schema.
        """
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        if config is None:
            config = lsst.meas.base.SingleFrameCompensatedGaussianFluxConfig()
        algorithm = lsst.meas.base.SingleFrameCompensatedGaussianFluxPlugin(
            config,
            "base_CompensatedGaussianFlux",
            schema,
            None
        )
        return algorithm, schema

    def _calcNumpyCompensatedGaussian(self, arr, var_arr, x_cent, y_cent, width, t):
        """Calculate the compensated Gaussian using numpy (for testing).

        Parameters
        ----------
        arr : `np.ndarray`
            Array of pixel values.
        var_arr : `np.ndarray`
            Array of variance values.
        x_cent : `float`
            x value of centroid.
        y_cent : `float`
            y value of centroid.
        width : `float`
            Width of inner kernel.
        t : `float`
            Scaling factor for outer kernel (outer_width = width*t).

        Returns
        -------
        flux : `float`
        variance : `float`
        """
        xx, yy = np.meshgrid(np.arange(arr.shape[0]), np.arange(arr.shape[1]))
        # Compute the inner and outer normalized Gaussian weights.
        inner = (1./(2.*np.pi*width**2.))*np.exp(-0.5*(((xx - x_cent)/width)**2. + ((yy - y_cent)/width)**2.))
        outer = (1./(2.*np.pi*(t*width)**2.))*np.exp(-0.5*(((xx - x_cent)/(t*width))**2.
                                                           + ((yy - y_cent)/(t*width))**2.))
        weight = inner - outer

        # Compute the weighted sum of the pixels.
        flux = np.sum(weight*arr)
        # And the normalization term, derived in Lupton et al. (in prep).
        flux *= 4.*np.pi*(width**2.)*(t**2. + 1)/(t**2. - 1)

        # And compute the variance
        variance = np.sum(weight*weight*var_arr)
        variance /= np.sum(weight*weight)

        # The variance normalization term, derived in Lupton et al. (in prep).
        variance *= 4.*np.pi*(width**2.)*(t**2 + 1)/t**2.

        return flux, variance

    def testCompensatedGaussianInnerProduct(self):
        """Test using the inner product routine directly."""

        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Extent2I(101, 101))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        psf_size = 2.0
        dataset.psfShape = lsst.afw.geom.Quadrupole(psf_size**2., psf_size**2., 0.0)
        centroid = lsst.geom.Point2D(50.0, 50.0)
        true_flux = 100000.0
        dataset.addSource(true_flux, centroid)

        # We need to set up the task in order to create the test dataset.
        task = self.makeSingleFrameMeasurementTask("base_CompensatedGaussianFlux")

        for width in self.all_widths:
            for t in self.all_ts:
                flux_0 = None
                var_0 = None
                for dc_offset in (0.0, -1.0, 1.0):
                    exposure, catalog = dataset.realize(10.0, task.schema, randomSeed=1000)

                    exposure.image.array += dc_offset

                    flux, var = _compensatedGaussianFiltInnerProduct(
                        exposure.image.array,
                        exposure.variance.array,
                        centroid.getX(),
                        centroid.getY(),
                        width,
                        t,
                    )

                    flux_numpy, var_numpy = self._calcNumpyCompensatedGaussian(
                        exposure.image.array,
                        exposure.variance.array,
                        centroid.getX(),
                        centroid.getY(),
                        width,
                        t,
                    )

                    # Compare values from c++ code to numpy code.
                    self.assertFloatsAlmostEqual(flux, flux_numpy, rtol=1e-6)
                    self.assertFloatsAlmostEqual(var, var_numpy, rtol=1e-6)

                    # If the kernel width is equal to the simulated PSF then it
                    # should be nearly unbiased.
                    if width == psf_size:
                        self.assertFloatsAlmostEqual(flux, true_flux, rtol=1e-3)

                    # And check biases with non-zero DC offset; these should be
                    # equal with some floating point tolerance.
                    if dc_offset == 0.0:
                        flux_0 = flux
                        var_0 = var
                    else:
                        self.assertFloatsAlmostEqual(flux, flux_0, rtol=1e-7)
                        self.assertFloatsAlmostEqual(var, var_0, rtol=1e-7)

    def testCompensatedGaussianSubPixels(self):
        """Test for correct instFlux as a function of sub-pixel position."""
        np.random.seed(12345)

        n_points = 100

        x_sub = np.random.uniform(low=0.0, high=1.0, size=n_points)
        y_sub = np.random.uniform(low=0.0, high=1.0, size=n_points)

        # We need to set up the task in order to create the test dataset.
        task = self.makeSingleFrameMeasurementTask("base_CompensatedGaussianFlux")

        for width in self.all_widths:
            for t in self.all_ts:
                fluxes = np.zeros(n_points)

                for i in range(n_points):
                    dataset = lsst.meas.base.tests.TestDataset(self.bbox_single)
                    centroid = lsst.geom.Point2D(50.0 + x_sub[i], 50.0 + y_sub[i])
                    dataset.addSource(50000.0, centroid)

                    exposure, catalog = dataset.realize(10.0, task.schema, randomSeed=i)
                    flux, var = _compensatedGaussianFiltInnerProduct(
                        exposure.image.array,
                        exposure.variance.array,
                        centroid.getX(),
                        centroid.getY(),
                        width,
                        t,
                    )

                    fluxes[i] = flux

                # Check for no correlation with x_sub and y_sub.
                fit_x = np.polyfit(x_sub, fluxes/50000.0, 1)
                self.assertLess(np.abs(fit_x[0]), 1e-3)
                fit_y = np.polyfit(y_sub, fluxes/50000.0, 1)
                self.assertLess(np.abs(fit_y[0]), 1e-3)

    def testCompensatedGaussianPlugin(self):
        """Test for correct instFlux given known position and shape.
        """
        # In the z-band, HSC images have a noise of about 40.0 ADU, and a background
        # offset of ~ -0.6 ADU/pixel.  This determines our test levels.
        for width in self.all_widths:
            for t in self.all_ts:
                for dc_offset in (0.0, -1.0, 1.0):
                    config = self.makeSingleFrameMeasurementConfig("base_CompensatedGaussianFlux")
                    config.algorithms["base_CompensatedGaussianFlux"].kernel_widths = [width]
                    config.algorithms["base_CompensatedGaussianFlux"].t = t

                    task = self.makeSingleFrameMeasurementTask(config=config)
                    exposure, catalog = self.dataset.realize(40.0, task.schema, randomSeed=0)
                    exposure.image.array += dc_offset
                    task.run(catalog, exposure)

                    filter_flux = catalog[f"base_CompensatedGaussianFlux_{width}_instFlux"]
                    filter_err = catalog[f"base_CompensatedGaussianFlux_{width}_instFluxErr"]
                    truth_flux = catalog["truth_instFlux"]

                    if width == self.psf_size:
                        # When the filter matches the PSF, we should get close to the true flux.
                        tol = np.sqrt((filter_err/filter_flux)**2. + 0.02**2.)
                        self.assertFloatsAlmostEqual(filter_flux, truth_flux, rtol=tol)
                    elif width > self.psf_size:
                        # When the filter is larger than the PSF, the filter flux will be
                        # greater than the truth flux.
                        np.testing.assert_array_less(truth_flux, filter_flux)

                    if dc_offset == 0.0:
                        # Use the no-offset run as a comparison for offset runs.
                        flux_0 = filter_flux
                    else:
                        # Note: this tolerance is determined empirically, but this is
                        # larger than preferable.
                        self.assertFloatsAlmostEqual(filter_flux, flux_0, rtol=5e-3)

                    # The ratio of the filter flux to the truth flux should be consistent.
                    # I'm not sure how to scale this with the error, so this is a loose
                    # tolerance now.
                    ratio = filter_flux / truth_flux
                    self.assertLess(np.std(ratio), 0.04)

    def testMonteCarlo(self):
        """Test an ideal simulation, with no noise.

        Demonstrate that:

        - We get exactly the right answer, and
        - The reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        nSamples = 500

        for width in self.larger_widths:
            for t in self.all_ts:
                config = lsst.meas.base.SingleFrameCompensatedGaussianFluxConfig()
                config.kernel_widths = [width]
                config.t = t

                algorithm, schema = self.makeAlgorithm(config=config)

                # Make a noiseless catalog.
                exposure, catalog = self.dataset_single.realize(1E-8, schema, randomSeed=1)

                # Only use the high-flux source for the error tests.
                record = catalog[0]
                algorithm.measure(record, exposure)
                inst_flux = record[f"base_CompensatedGaussianFlux_{width}_instFlux"]

                for noise in (0.001, 0.01, 0.1):
                    fluxes = np.zeros(nSamples)
                    errs = np.zeros_like(fluxes)

                    for repeat in range(nSamples):
                        # By using ``repeat`` to seed the RNG, we get results which
                        # fall within the tolerances defined below. If we allow this
                        # test to be truly random, passing becomes RNG-dependent.
                        exposure_samp, catalog_samp = self.dataset_single.realize(
                            noise*inst_flux,
                            schema,
                            randomSeed=repeat,
                        )
                        record_samp = catalog_samp[0]
                        algorithm.measure(record_samp, exposure_samp)
                        fluxes[repeat] = record_samp[f"base_CompensatedGaussianFlux_{width}_instFlux"]
                        errs[repeat] = record_samp[f"base_CompensatedGaussianFlux_{width}_instFluxErr"]

                    err_mean = np.mean(errs)
                    flux_std = np.std(fluxes)
                    self.assertFloatsAlmostEqual(err_mean, flux_std, rtol=0.10)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
