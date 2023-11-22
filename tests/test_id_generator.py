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
from typing import Callable

from lsst.pex.config import Config
from lsst.afw.table import SourceTable
from lsst.daf.butler import DimensionUniverse, DataCoordinate
from lsst.meas.base import (
    IdGenerator,
    DetectorVisitIdGeneratorConfig,
    DetectorExposureIdGeneratorConfig,
    SkyMapIdGeneratorConfig,
)


class _TestConfig(Config):
    skymap = SkyMapIdGeneratorConfig.make_field()
    visit = DetectorVisitIdGeneratorConfig.make_field()
    exposure = DetectorExposureIdGeneratorConfig.make_field()


class IdGeneratorTestCase(unittest.TestCase):
    """Tests for the IdGenerator class and its config-based construction
    pattern.
    """

    def setUp(self):
        self.schema = SourceTable.makeMinimalSchema()
        self.universe = DimensionUniverse()

    def test_no_packer(self):
        """Test a simple IdGenerator that doesn't embed anything."""
        id_gen = IdGenerator()
        self.check_invariants(id_gen)

    def test_visit(self):
        """Test an IdGenerator that packs {visit, detector} data IDs."""
        data_id = DataCoordinate.standardize(
            instrument="I",
            visit=312,
            detector=5,
            universe=self.universe,
        )
        config = _TestConfig()
        config.visit.packer.name = "observation"
        config.visit.packer["observation"].n_observations = 10000
        config.visit.packer["observation"].n_detectors = 99
        config.visit.n_releases = 8
        config.visit.release_id = 2
        id_generator = config.visit.apply(data_id)
        self.check_invariants(
            id_generator,
            unpacker=IdGenerator.unpacker_from_config(
                config.visit, data_id.subset(self.universe.conform(["instrument"]))
            ),
            expected_release_id=2,
        )

    def test_exposure(self):
        """Test an IdGenerator that packs {exposure, detector} data IDs."""
        data_id = DataCoordinate.standardize(
            instrument="I",
            exposure=312,
            detector=5,
            universe=self.universe,
        )
        config = _TestConfig()
        config.exposure.packer.name = "observation"
        config.exposure.packer["observation"].n_observations = 10000
        config.exposure.packer["observation"].n_detectors = 99
        config.exposure.n_releases = 4
        config.exposure.release_id = 3
        id_generator = config.exposure.apply(data_id)
        self.check_invariants(
            id_generator,
            unpacker=IdGenerator.unpacker_from_config(
                config.exposure, data_id.subset(self.universe.conform(["instrument"]))
            ),
            expected_release_id=3,
        )

    def test_skymap(self):
        """Test an IdGenerator that packs {tract, patch, band} data IDs."""
        config = _TestConfig()
        config.skymap.packer.n_tracts = 11
        config.skymap.packer.n_patches = 9
        data_id = DataCoordinate.standardize(
            skymap="S", tract=9, patch=5, band="r", universe=self.universe
        )
        id_generator = config.skymap.apply(data_id)
        self.check_invariants(
            id_generator,
            unpacker=IdGenerator.unpacker_from_config(
                config.skymap, data_id.subset(self.universe.conform(["skymap"]))
            ),
        )

    def check_invariants(
        self,
        id_gen: IdGenerator,
        unpacker: Callable[[int], tuple[DataCoordinate, int]] | None = None,
        expected_release_id: int = 0,
    ):
        """Check methods of the `IdGenerator` class for self-consistency and
        expected values.
        """
        catalog = id_gen.make_source_catalog(self.schema)
        for _ in range(5):
            catalog.addNew()
        array = id_gen.arange(1, 6)
        self.assertEqual(list(catalog["id"]), list(array))
        if unpacker is not None:
            expected_data_id = id_gen.data_id
            for i in range(5):
                embedded_release_id, embedded_data_id, counter = unpacker(array[i])
                self.assertEqual(counter, i + 1)
                self.assertEqual(embedded_data_id, expected_data_id)
                self.assertEqual(embedded_release_id, expected_release_id)


if __name__ == "__main__":
    unittest.main()
