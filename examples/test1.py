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

import math
import sys
import os
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementTask
from lsst.meas.base.base import *
from lsst.meas.base.tests import *
import unittest
import lsst.utils.tests
import numpy

numpy.random.seed(1234)


if __name__ == "__main__":
    if not len(sys.argv) == 2:
        print "Usage: %s visit"%(sys.argv[0],)
        sys.exit(1)
    visit = sys.argv[1] 
    DATA_FILE = "mmout1/src/v%s-fi/R22/S11.fits"%(visit,)
    DATA_FILE0 = "mmout0/src/v%s-fi/R22/S11.fits"%(visit,)

    measCat0 = SourceCatalog.readFits(DATA_FILE0)
    measCat = SourceCatalog.readFits(DATA_FILE)
    print measCat.getCentroidDefinition()
    print measCat.getPsfFluxDefinition()
    print measCat.getModelFluxDefinition()
    print measCat0.getCentroidDefinition()
    print measCat0.getPsfFluxDefinition()
    print measCat0.getModelFluxDefinition()
    assert(len(measCat) == len(measCat0))
    for i in range(len(measCat)):
        record = measCat[i]
        record0 = measCat0[i]
        error = record.getCentroidErr()
        error0 = record0.getCentroidErr()
        flag = record.getCentroidFlag()
        flag0 = record0.getCentroidFlag()
        value = record.getCentroid()
        value0 = record0.getCentroid()
        label = "Centroid: "
        if not (value==value0) and not(flag and flag0):
            print label, record.getCentroid(), record.getId(), record.getCentroid(), record0.getCentroid(), record.getCentroidFlag(), record0.getCentroidFlag()

        value = record.getPsfFlux()
        value0 = record0.getPsfFlux()
        error = record.getPsfFluxErr()
        error0 = record0.getPsfFluxErr()
        flag = record.getPsfFluxFlag()
        flag0 = record0.getPsfFluxFlag()
        label = "PsfFlux: "
        if not (abs((value-value0)/value0)<.001) and not(numpy.isnan(value)) and not record.get("base_PsfFlux_flag_edge"):
            print label, "Values: ", record.getCentroid(), record.getId(), record.getPsfFlux(), record0.getPsfFlux(), record.getPsfFluxFlag(), record0.getPsfFluxFlag()
        if (flag0 != flag) and not record.get("base_PsfFlux_flag_edge"):
            print label, "Flags: ", record.getCentroid(), record.getId(), record.getPsfFlux(), record0.getPsfFlux(), record.getPsfFluxFlag(), record0.getPsfFluxFlag()
        if not (abs((error-error0)/error0)<.1) and not(numpy.isnan(error)) and not record.get("base_PsfFlux_flag_edge"):
            print label, "Errors: ", record.getCentroid(), record.getId(), record.getPsfFluxErr(), record0.getPsfFluxErr(), record.getPsfFluxFlag(), record0.getPsfFluxFlag()

        value = record.getInstFlux()
        value0 = record0.getInstFlux()
        error = record.getInstFluxErr()
        error0 = record0.getInstFluxErr()
        flag = record.getInstFluxFlag()
        flag0 = record0.getInstFluxFlag()
        error = record.getInstFluxErr()
        error0 = record0.getInstFluxErr()
        label = "NaiveFlux: "
        if not (value==value0) and not(numpy.isnan(value) and numpy.isnan(value0)):
            print label, record.getCentroid(), record.getId(), record.getInstFlux(), record0.getInstFlux(), record.getInstFluxFlag(), record0.getInstFluxFlag()
        if (flag0 != flag):
            print label, "Flags: ", record.getCentroid(), record.getId(), record.getInstFlux(), record0.getInstFlux(), record.getInstFluxFlag(), record0.getInstFluxFlag()
        if not (abs((error-error0)/error0)<.0001) and not(numpy.isnan(error)):
            print label, "Errors: ", record.getCentroid(), record.getId(), record.getInstFluxErr(), record0.getInstFluxErr(), record.getInstFluxFlag(), record0.getInstFluxFlag()

        value = record.getApFlux()
        value0 = record0.getApFlux()
        error = record.getApFluxErr()
        error0 = record0.getApFluxErr()
        flag = record.getApFluxFlag()
        flag0 = record0.getApFluxFlag()
        error = record.getApFluxErr()
        error0 = record0.getApFluxErr()
        label = "SincFlux: "
        if not (value==value0) and not(numpy.isnan(value) and numpy.isnan(value0)):
            print label, record.getCentroid(), record.getId(), record.getApFlux(), record0.getApFlux(), record.getApFluxFlag(), record0.getApFluxFlag()
        if (flag0 != flag):
            print label, "Flags: ", record.getCentroid(), record.getId(), record.getApFlux(), record0.getApFlux(), record.getApFluxFlag(), record0.getApFluxFlag()
        if not (abs((error-error0)/error0)<.0001) and not(numpy.isnan(error)):
            print label, "Errors: ", record.getCentroid(), record.getId(), record.getApFluxErr(), record0.getApFluxErr(), record.getApFluxFlag(), record0.getApFluxFlag()

        value = record.getModelFlux()
        value0 = record0.getModelFlux()
        error = record.getModelFluxErr()
        error0 = record0.getModelFluxErr()
        flag = record.getModelFluxFlag()
        flag0 = record0.getModelFluxFlag()
        flag0 = record0.getModelFluxFlag()
        error = record.getModelFluxErr()
        label = "GaussianFlux: "
        if not (value==value0) and not(numpy.isnan(value) and numpy.isnan(value0)):
            print label, record.getCentroid(), record.getId(), record.getModelFlux(), record0.getModelFlux(), record.getModelFluxFlag(), record0.getModelFluxFlag()
        if (flag0 != flag):
            print label, "Flags: ", record.getCentroid(), record.getId(), record.getModelFlux(), record0.getModelFlux(), record.getModelFluxFlag(), record0.getModelFluxFlag()
        if not (abs((error-error0)/error0)<.0001) and not(numpy.isnan(error)):
            print label, "Errors: ", record.getCentroid(), record.getId(), record.getModelFluxErr(), record0.getModelFluxErr(), record.getModelFluxFlag(), record0.getModelFluxFlag()

