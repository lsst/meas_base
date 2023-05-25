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

import lsst.meas.base as measBase
import lsst.meas.base.tests
import lsst.utils.tests
import numpy as np


class MomentsClassificationTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -20),
                                    lsst.geom.Extent2I(250, 150))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)

        self.n_stars = 1
        self.n_gals = 1
        # First 10 sources are point sources
        self.dataset.addSource(1000.0, lsst.geom.Point2D(50.1, 49.8))
        # Following 10 sources are extended sources
        self.dataset.addSource(5000.0, lsst.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.Quadrupole(1, 1.2, 0.3))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def testSingleFramePlugin(self):
        config = measBase.SingleFrameMeasurementConfig()
        task = self.makeSingleFrameMeasurementTask(config=config)
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=3)
        task.run(catalog, exposure)
        for ii in range(self.n_stars):
            self.assertLess(catalog[ii].get("base_ClassificationSizeExtendedness_value"), 0.1)
        for ii in range(self.n_stars, self.n_stars + self.n_gals):
            self.assertGreater(catalog[ii].get("base_ClassificationSizeExtendedness_value"), 0.02)

    @lsst.utils.tests.methodParameters(noise=(0.001, 0.01))
    def testMonteCarlo(self, noise: float, n_trials: int = 100):
        """Test an ideal simulation, with no noise.

        Demonstrate that:

        - We get exactly the right answer, and
        - The reported uncertainty agrees with a Monte Carlo test of the noise.

        Parameters
        ----------
        noise : float
            Noise level to use in the simulation.
        n_trials : int
            Number of trials to use in the Monte Carlo test.
        """
        config = measBase.SingleFrameMeasurementConfig()
        task = self.makeSingleFrameMeasurementTask(config=config)

        star_measures, galaxy_measures = [], []
        for ii in range(n_trials):
            exposure, catalog = self.dataset.realize(1000.0*noise, task.schema, randomSeed=ii)
            task.run(catalog, exposure)
            for ii in range(self.n_stars):
                star_measures.append(catalog[ii].get("base_ClassificationSizeExtendedness_value"))
            for ii in range(self.n_stars, self.n_stars + self.n_gals):
                galaxy_measures.append(catalog[ii].get("base_ClassificationSizeExtendedness_value"))

        # Mapping noise level to thresholds for stars and galaxies
        star_threshold = {
            0.001: 0.01,
            0.01: 0.1
        }
        galaxy_threshold = {
            0.001: 0.25,
            0.01: 0.25,
        }
        self.assertLess(np.mean(star_measures), star_threshold[noise])
        self.assertLess(np.percentile(star_measures, 50), star_threshold[noise])
        self.assertGreater(np.mean(galaxy_measures), galaxy_threshold[noise])
        self.assertGreater(np.percentile(galaxy_measures, 50), galaxy_threshold[noise])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
