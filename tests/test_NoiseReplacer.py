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
import lsst.afw.detection
import lsst.afw.table
import lsst.meas.base.tests
import lsst.utils.tests
from lsst.daf.base import PropertyList
from lsst.utils.tests import methodParameters


@lsst.meas.base.register("test_NoiseReplacer")
class NoiseReplacerTestPlugin(lsst.meas.base.SingleFramePlugin):
    """Plugin that sums instFlux inside and outside the source footprint.
    """

    @staticmethod
    def getExecutionOrder():
        return 2.0

    def __init__(self, config, name, schema, metadata):
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.insideKey = schema.addField("%s_inside" % (name,), type=np.float64,
                                         doc="instFlux inside footprint")
        self.outsideKey = schema.addField("%s_outside" % (name,), type=np.float64,
                                          doc="instFlux outside footprint")

    def measure(self, measRecord, exposure):
        footprint = measRecord.getFootprint()
        fullArray = exposure.image.array
        insideArray = np.zeros(footprint.getArea(), dtype=fullArray.dtype)
        footprint.spans.flatten(insideArray, fullArray, exposure.getXY0())
        insideFlux = float(insideArray.sum())
        outsideFlux = float(fullArray.sum()) - insideFlux
        measRecord.set(self.insideKey, insideFlux)
        measRecord.set(self.outsideKey, outsideFlux)


class NoiseReplacerTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(240, 260))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(120000.0, lsst.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.Quadrupole(8, 9, 3))
        with self.dataset.addBlend() as family:
            family.addChild(110000.0, lsst.geom.Point2D(65.2, 150.7),
                            lsst.afw.geom.Quadrupole(7, 5, -1))
            family.addChild(140000.0, lsst.geom.Point2D(72.3, 149.1))
            family.addChild(90000.0, lsst.geom.Point2D(68.5, 156.9))

    @methodParameters(noiseSource=['measure', 'variance', 'meta', 'variance_median'],
                      variance=[1.0, 1.0, 2.0, 1.0])
    def testNoiseReplacer(self, noiseSource, variance):
        """Test noise replacement in SFM with ''measure'' mode.

        We compare the instFlux inside and outside source Footprints on an
        extremely high S/N image.
        """
        # We choose a random seed which causes the test to pass.
        task = self.makeSingleFrameMeasurementTask("test_NoiseReplacer")
        task.config.noiseReplacer.noiseSource = noiseSource
        exposure, catalog = self.dataset.realize(variance, task.schema, randomSeed=0)
        if noiseSource == 'meta':
            md = PropertyList()
            md['BGMEAN'] = variance
            exposure.setMetadata(md)
        task.run(catalog, exposure)
        sumVariance = exposure.variance.array.sum()
        for record in catalog:
            self.assertFloatsAlmostEqual(record.get("test_NoiseReplacer_inside"),
                                         record.get("truth_instFlux"), rtol=1E-3)

            # N.B. Next line checks that a random value is correct to a
            # statistical 1-sigma prediction; some RNG seeds may cause it to
            # fail (indeed, 67% should)
            self.assertLess(record.get("test_NoiseReplacer_outside"), np.sqrt(sumVariance))

    def tearDown(self):
        del self.bbox
        del self.dataset


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
