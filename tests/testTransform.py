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

import unittest

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.base as measBase
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExcept
import lsst.utils.tests
import testLib

try:
    type(verbose)
except NameError:
    verbose = 0


def makeWcs():
    """Provide a simple WCS for use in testing"""
    # The parameters given are placeholders; their values are unimportant
    md = dafBase.PropertySet()
    md.set("NAXIS", 2)
    md.set("CTYPE1", "RA---TAN")
    md.set("CTYPE2", "DEC--TAN")
    md.set("CRPIX1", 0)
    md.set("CRPIX2", 0)
    md.set("CRVAL1", 0)
    md.set("CRVAL2", 0)
    md.set("RADECSYS", "FK5")
    md.set("EQUINOX", 2000.0)
    return afwImage.makeWcs(md)


@pexConfig.wrap(testLib.SillyCentroidControl)
class SillyCentroidConfig(pexConfig.Config):
    pass


class TransformTestCase(lsst.utils.tests.TestCase):
    pluginName = "base_SillyCentroid"
    centroidPosition = (1.0, -1.0)

    def _generateCatalog(self):
        """Returns a SourceCatalog with one entry"""
        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField(self.pluginName + "_x", type=float)
        schema.addField(self.pluginName + "_y", type=float)
        cat = afwTable.SourceCatalog(schema)
        source = cat.addNew()
        source.set(self.pluginName + "_x", self.centroidPosition[0])
        source.set(self.pluginName + "_y", self.centroidPosition[1])
        return cat

    def _performTransform(self, transformClass, inCat, doExtend=True):
        """Operate on inCat with a transform of class transformClass"""
        mapper = afwTable.SchemaMapper(inCat.schema)
        config = SillyCentroidConfig()
        transform = transformClass(config, self.pluginName, mapper)
        outCat = afwTable.BaseCatalog(mapper.getOutputSchema())
        if doExtend:
            outCat.extend(inCat, mapper=mapper)
        transform(inCat, outCat, makeWcs(), afwImage.Calib())
        return outCat

    def _checkSillyOutputs(self, inCat, outCat):
        """Check that outCat corresponds to inCat after application of SillyTransform"""
        # SillyTransform looks for fields named ``name_x`` and ``name_y``; it
        # copies them to the output, and adds ``name_reverse_x`` and ``_y``
        # which are equal to the corresponding inputs multiplied by -1.
        for inSrc, outSrc in zip(inCat, outCat):
            # The source x, y should have been copied to the output table
            self.assertEqual(outSrc[self.pluginName + "_x"], inSrc[self.pluginName + "_x"])
            self.assertEqual(outSrc[self.pluginName + "_y"], inSrc[self.pluginName + "_y"])

            # And the reversed position added
            self.assertEqual(outSrc[self.pluginName + "_reverse_x"], -1.0 * inSrc[self.pluginName + "_x"])
            self.assertEqual(outSrc[self.pluginName + "_reverse_y"], -1.0 * inSrc[self.pluginName + "_y"])

        # Other entries from the source table have not been copied
        for name in ("id", "coord_ra", "coord_dec", "parent"):
            self.assertIn(name, inCat.schema)
            self.assertNotIn(name, outCat.schema)

    def testNullTransform(self):
        """The NullTransform passes through nothing"""
        inCat = self._generateCatalog()
        with self.assertRaises(pexExcept.LengthError):
            self._performTransform(measBase.NullTransform, inCat, False)
        outCat = self._performTransform(measBase.NullTransform, inCat)
        self.assertEqual(len(inCat), len(outCat))
        self.assertEqual(outCat.schema.getFieldCount(), 0)

    def testPassThroughTransform(self):
        """The PassThroughTransform copies all fields starting with the plugin name"""
        inCat = self._generateCatalog()
        with self.assertRaises(pexExcept.LengthError):
            self._performTransform(measBase.PassThroughTransform, inCat, False)
        outCat = self._performTransform(measBase.PassThroughTransform, inCat)
        self.assertEqual(len(inCat), len(outCat))
        for inSrc, outSrc in zip(inCat, outCat):
            for fieldname in inCat.schema.extract(self.pluginName + "*").keys():
                self.assertEqual(inSrc.get(fieldname), outSrc.get(fieldname))

    def testPythonConfig(self):
        """The Python Config should be automatically converted to Control when calling a C++ transform."""
        inCat = self._generateCatalog()
        with self.assertRaises(pexExcept.LengthError):
            self._performTransform(testLib.SillyTransform, inCat, False)
        outCat = self._performTransform(testLib.SillyTransform, inCat)
        self._checkSillyOutputs(inCat, outCat)

    def testApplyCppTransform(self):
        """Test that we can apply a simple C++ transform"""
        inCat = self._generateCatalog()
        sillyControl = testLib.SillyCentroidControl()
        mapper = afwTable.SchemaMapper(inCat.schema)
        sillyTransform = testLib.SillyTransform(sillyControl, self.pluginName, mapper)
        outCat = afwTable.BaseCatalog(mapper.getOutputSchema())
        outCat.extend(inCat, mapper=mapper)
        self.assertEqual(len(inCat), len(outCat))
        sillyTransform(inCat, outCat, makeWcs(), afwImage.Calib())
        self._checkSillyOutputs(inCat, outCat)


class AlgorithmConfigurationTestCase(lsst.utils.tests.TestCase):

    def testDefaultTransform(self):
        """By default, we perform no transformations"""
        self.assertEqual(measBase.BasePlugin.getTransformClass(), measBase.PassThroughTransform)

    def testWrapAlgorithm(self):
        """Test that the appropriate transform is provided for wrapped algorithms"""
        # By default, we inherit from BasePlugin
        # NB the choice of algorithm and executionOrder is arbitrary
        singleFrame, forced = measBase.wrapSimpleAlgorithm(measBase.PsfFluxAlgorithm,
                                                           Control=measBase.PsfFluxControl,
                                                           executionOrder=measBase.BasePlugin.FLUX_ORDER,
                                                           doRegister=False)
        self.assertEqual(singleFrame.getTransformClass(), measBase.BasePlugin.getTransformClass())
        # Unless we override
        singleFrame, forced = measBase.wrapSimpleAlgorithm(measBase.PsfFluxAlgorithm,
                                                           Control=measBase.PsfFluxControl,
                                                           executionOrder=measBase.BasePlugin.FLUX_ORDER,
                                                           doRegister=False,
                                                           TransformClass=measBase.PassThroughTransform)
        self.assertEqual(singleFrame.getTransformClass(), measBase.PassThroughTransform)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
