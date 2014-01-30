#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import os
import unittest
import numpy

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.base

class ResultMappersTestCase(lsst.utils.tests.TestCase):
    """Test case for utility classes used in plugin wrappers
    """

    def testFluxResultMapper(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        mapper = lsst.meas.base.FluxResultMapper(schema, "test", lsst.meas.base.DIAGONAL_ONLY)
        table = lsst.afw.table.SourceTable.make(schema)
        result = lsst.meas.base.FluxResult()
        result.flux = 34.5
        result.fluxSigma = 3.25
        record = table.makeRecord()
        mapper.apply(record, result)
        self.assertEqual(record.get("test_flux"), result.flux)
        self.assertEqual(record.get("test_fluxSigma"), result.fluxSigma)
        self.assertEqual(record.get("test_flag"), False)
        mapper.fail(record)
        self.assertEqual(record.get("test_flag"), True)

    def testCentroidResultMapper(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        mapper = lsst.meas.base.CentroidResultMapper(schema, "test", lsst.meas.base.DIAGONAL_ONLY)
        table = lsst.afw.table.SourceTable.make(schema)
        result = lsst.meas.base.CentroidResult()
        result.x = 1.2
        result.y = 3.4
        result.xSigma = 2.0
        result.ySigma = 3.0
        record = table.makeRecord()
        mapper.apply(record, result)
        self.assertEqual(record.get("test_x"), result.x)
        self.assertEqual(record.get("test_y"), result.y)
        self.assertEqual(record.get("test_xSigma"), result.xSigma)
        self.assertEqual(record.get("test_ySigma"), result.ySigma)
        self.assertEqual(record.get("test_flag"), False)
        mapper.fail(record)
        self.assertEqual(record.get("test_flag"), True)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ResultMappersTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
