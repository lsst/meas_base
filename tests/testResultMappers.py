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

    def testFluxAlgorithmMapper(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        mapper = lsst.meas.base.FluxAlgorithmMapper(schema, "test")
        table = lsst.afw.table.SourceTable.make(schema)
        result = lsst.meas.base.FluxAlgorithmResult(34.5, 3.25)
        record = table.makeRecord()
        mapper.apply(record, result)
        self.assertEqual(record.get("test.value"), result.value)
        self.assertEqual(record.get("test.err"), result.err)
        self.assertEqual(record.get("test.flag"), False)
        mapper.fail(record)
        self.assertEqual(record.get("test.flag"), True)

    def testCentroidAlgorithmMapper(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        mapper = lsst.meas.base.CentroidAlgorithmMapper(schema, "test")
        table = lsst.afw.table.SourceTable.make(schema)
        result = lsst.meas.base.CentroidAlgorithmResult(lsst.afw.geom.Point2D(1.2, 3.4),
                                                        numpy.array([[2.0, 0.1], [0.1, 3.0]],
                                                                    dtype=numpy.float32))
        record = table.makeRecord()
        mapper.apply(record, result)
        self.assertEqual(record.get("test.value"), result.value)
        self.assertClose(record.get("test.cov"), result.cov, rtol=0.0, atol=0.0)
        self.assertEqual(record.get("test.flag"), False)
        mapper.fail(record)
        self.assertEqual(record.get("test.flag"), True)


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
