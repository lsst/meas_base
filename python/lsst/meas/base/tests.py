#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2014 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import lsst.utils.tests

import lsst.afw.table
import lsst.afw.image

from .sfm import *
from .forcedImage import *

DATA_DIR = os.path.join(os.path.split(__file__)[0], "../../../../tests/data")

class AlgorithmTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.truth = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-01.fits"))
        self.calexp = lsst.afw.image.ExposureF.readFits(os.path.join(DATA_DIR, "calexp-01.fits"))
        self.fluxKey = self.truth.schema.find("truth.flux").key
        self.centroidKey = self.truth.schema.find("truth.centroid").key
        self.shapeKey = self.truth.schema.find("truth.shape").key

    def tearDown(self):
        del self.truth
        del self.calexp
        del self.fluxKey
        del self.centroidKey

    def runSingleFrameMeasurementTask(self, plugin, dependencies=(), flags=None, config=None):
        if config is None:
            config = SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = None
        config.slots.shape = None
        config.slots.modelFlux = None
        config.slots.apFlux = None
        config.slots.psfFlux = None
        config.slots.instFlux = None
        config.plugins.names = (plugin,) + tuple(dependencies)
        schemaMapper = lsst.afw.table.SchemaMapper(self.truth.schema)
        schemaMapper.addMinimalSchema(self.truth.schema)
        task = SingleFrameMeasurementTask(schema=schemaMapper.editOutputSchema(), config=config, flags=flags)
        measCat = lsst.afw.table.SourceCatalog(task.schema)
        measCat.extend(self.truth, schemaMapper)
        measCat.getTable().defineModelFlux(self.fluxKey)
        measCat.getTable().defineCentroid(self.centroidKey)
        measCat.getTable().defineShape(self.shapeKey)
        task.run(measCat, self.calexp)
        return measCat
