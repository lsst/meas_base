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

import unittest

import lsst.utils.tests
import lsst.meas.base
import lsst.meas.base.tests
import lsst.afw.table
from lsst.meas.base.baseLib import MeasurementError
from lsst.meas.base import FlagDefinition, FlagDefinitionVector, FlagHandler
from lsst.meas.base.tests import (AlgorithmTestCase)

class FlagHandlerTestCase(AlgorithmTestCase):
    
    def testFlagHandler(self):
        control = lsst.meas.base.GaussianCentroidControl()
        alg = lsst.meas.base.GaussianCentroidAlgorithm
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        plugin = alg(control, 'test', schema)
        cat = lsst.afw.table.SourceCatalog(schema)
        schema = cat.getSchema()
        subSchema = schema["test"]

        # This is a FlagDefinition structure like a plugin might have
        FAILURE = 0
        FIRST = 1
        SECOND = 2
        flagDefs = [ FlagDefinition("General Failure", "general failure error"),
            FlagDefinition("1st error", "this is the first failure type"),
            FlagDefinition("2nd error", "this is the second failure type")
        ]
        fh = FlagHandler.addFields(schema, "test",
            FlagDefinitionVector(flagDefs))

        # Check to be sure that the FlagHandler was correctly initialized
        for index, flagDef in enumerate(flagDefs):
           assert(flagDef.name == flagDefs[index].name)
           assert(flagDef.doc == flagDefs[index].doc)

        catalog = lsst.afw.table.SourceCatalog(schema)

        # Now check to be sure that all of the known failures set the bits correctly
        record = catalog.addNew()
        fh.handleFailure(record)
        self.assertTrue(fh.getValue(record, 0))
        self.assertFalse(fh.getValue(record, 1))
        self.assertFalse(fh.getValue(record, 2))
        record = catalog.addNew()

        error = MeasurementError(fh.getDefinition(FAILURE).doc, FAILURE)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, 0))
        self.assertFalse(fh.getValue(record, 1))
        self.assertFalse(fh.getValue(record, 2))

        record = catalog.addNew()
        error = MeasurementError(fh.getDefinition(FIRST).doc, FIRST)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, 0))
        self.assertTrue(fh.getValue(record, 1))
        self.assertFalse(fh.getValue(record, 2))

        record = catalog.addNew()
        error = MeasurementError(fh.getDefinition(SECOND).doc, SECOND)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, 0))
        self.assertFalse(fh.getValue(record, 1))
        self.assertTrue(fh.getValue(record, 2))

        record = catalog.addNew()
        error = MeasurementError("Custom error message", FIRST)
        fh.handleFailure(record, error.cpp)
        self.assertTrue(fh.getValue(record, 0))
        self.assertTrue(fh.getValue(record, 1))
        self.assertFalse(fh.getValue(record, 2))

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(FlagHandlerTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
