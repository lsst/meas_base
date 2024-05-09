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


class CompensatedTophatFluxTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
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

        self.all_apertures = (12, 14)
        self.all_scales = ((1.2, 1.7), (1.5, 2.0))

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
            config = lsst.meas.base.SingleFrameCompensatedTophatFluxConfig()
        algorithm = lsst.meas.base.SingleFrameCompensatedTophatFluxPlugin(
            config,
            "base_CompensatedTophatFlux",
            schema,
            None
        )
        return algorithm, schema

    def testCompensatedTophatPlugin(self):
        """Test for correct instFlux given known position and shape.
        """
        # In the z-band, HSC images have a noise of about 40.0 ADU, and a background
        # offset of ~ -0.6 ADU/pixel.  This determines our test levels.
        for aperture in self.all_apertures:
            for inner_scale, outer_scale in self.all_scales:
                for dc_offset in (0.0, -1.0, 1.0):
                    config = self.makeSingleFrameMeasurementConfig("base_CompensatedTophatFlux")
                    config.algorithms["base_CompensatedTophatFlux"].apertures = [aperture]
                    config.algorithms["base_CompensatedTophatFlux"].inner_scale = inner_scale
                    config.algorithms["base_CompensatedTophatFlux"].outer_scale = outer_scale

                    task = self.makeSingleFrameMeasurementTask(config=config)
                    exposure, catalog = self.dataset.realize(40.0, task.schema, randomSeed=0)
                    exposure.image.array += dc_offset
                    task.run(catalog, exposure)

                    filter_flux = catalog[f"base_CompensatedTophatFlux_{aperture}_instFlux"]
                    filter_err = catalog[f"base_CompensatedTophatFlux_{aperture}_instFluxErr"]
                    truth_flux = catalog["truth_instFlux"]

                    np.testing.assert_array_less(truth_flux, filter_flux + 3*filter_err)
                    np.testing.assert_array_less(filter_flux - 3*filter_err, truth_flux)

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

    def testCompensatedTophatPluginFailure(self):
        """Test that the correct flag is set on failures."""
        config = self.makeSingleFrameMeasurementConfig("base_CompensatedTophatFlux")
        config.algorithms["base_CompensatedTophatFlux"].apertures = [5, 15]

        task = self.makeSingleFrameMeasurementTask(config=config)
        exposure, catalog = self.dataset.realize(40.0, task.schema, randomSeed=0)
        # Modify two objects to trigger the 2 failure conditions.
        catalog[0]["slot_Centroid_x"] = -20.0
        catalog[1]["slot_Centroid_x"] = -10.5
        task.run(catalog, exposure)
        self.assertTrue(catalog["base_CompensatedTophatFlux_5_flag"][0])
        self.assertTrue(catalog["base_CompensatedTophatFlux_5_flag_bounds"][0])
        self.assertTrue(catalog["base_CompensatedTophatFlux_15_flag"][0])
        self.assertTrue(catalog["base_CompensatedTophatFlux_15_flag_bounds"][0])
        self.assertTrue(catalog["base_CompensatedTophatFlux_5_flag"][1])
        self.assertTrue(catalog["base_CompensatedTophatFlux_5_flag_bounds"][1])

    def testMonteCarlo(self):
        """Test an ideal simulation, with no noise.

        Demonstrate that:

        - We get exactly the right answer, and
        - The reported uncertainty agrees with a Monte Carlo test of the noise.
        """
        nSamples = 500

        for aperture in self.all_apertures:
            for inner_scale, outer_scale in self.all_scales:
                config = lsst.meas.base.SingleFrameCompensatedTophatFluxConfig()
                config.apertures = [aperture]
                config.inner_scale = inner_scale
                config.outer_scale = outer_scale

                algorithm, schema = self.makeAlgorithm(config=config)

                # Make a noiseless catalog.
                exposure, catalog = self.dataset_single.realize(1E-8, schema, randomSeed=1)

                # Only use the high-flux source for the error tests.
                record = catalog[0]
                algorithm.measure(record, exposure)
                inst_flux = record[f"base_CompensatedTophatFlux_{aperture}_instFlux"]

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
                        fluxes[repeat] = record_samp[f"base_CompensatedTophatFlux_{aperture}_instFlux"]
                        errs[repeat] = record_samp[f"base_CompensatedTophatFlux_{aperture}_instFluxErr"]

                    err_mean = np.mean(errs)
                    flux_std = np.std(fluxes)
                    self.assertFloatsAlmostEqual(err_mean, flux_std, rtol=0.03)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
