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
import os
from lsst.afw.table import Schema,SchemaMapper,SourceCatalog,SourceTable
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin, SingleFrameMeasurementTask
from lsst.meas.base.base import *
from lsst.meas.base.tests import *
import unittest
import lsst.utils.tests
import numpy

numpy.random.seed(1234)


DATA_FILE = os.path.join(os.environ["work"], "mmout12/src/v100-fi/R22/S11.fits")
DATA_FILE0 = os.path.join(os.environ["work"], "mmout11/src/v100-fi/R22/S11.fits")

if __name__ == "__main__":
#  Read a catalog for the old (meas_algorithm) and new (meas_base) algorithms
#  These are the result of a complete run of processCcd with default configurations of measurement.py
#  and sfm.py.  "0" means the old meas_algorithm algorithms.  Only the measurement task is different.
#  The pipeline preparatory to measurement should be all the same for both catalogs.   

    errorLimit = .01
    valueLimit = .001
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
        # Check the Centroids, but only if at least one of the flags is false.  
        # The centroids behave differently on error for the old and new algorithms 
        if not (value==value0) and not(flag and flag0):
            print label, record.getCentroid(), record.getId(), record.getCentroid(), record0.getCentroid(), record.getCentroidFlag(), record0.getCentroidFlag()

        value = record.getPsfFlux()
        value0 = record0.getPsfFlux()
        error = record.getPsfFluxErr()
        error0 = record0.getPsfFluxErr()
        flag = record.getPsfFluxFlag()
        flag0 = record0.getPsfFluxFlag()
        label = "PsfFlux: "
        if not (abs((value-value0)/value0)<valueLimit) and not(numpy.isnan(value)) and not record.get("base_PsfFlux_flag_edge"):
            print label, "Values differ: ", record.getCentroid(), record.getId(), record.getPsfFlux(), record0.getPsfFlux(), record.getPsfFluxFlag(), record0.getPsfFluxFlag()
        if (flag0 != flag) and not record.get("base_PsfFlux_flag_edge"):
            print label, "Flags differ: ", record.getCentroid(), record.getId(), record.getPsfFlux(), record0.getPsfFlux(), record.getPsfFluxFlag(), record0.getPsfFluxFlag()
        if not (abs((error-error0)/error0)<errorLimit) and not(numpy.isnan(error)) and not record.get("base_PsfFlux_flag_edge"):
            print label, "Errors differ: ", record.getCentroid(), record.getId(), record.getPsfFluxErr(), record0.getPsfFluxErr(), record.getPsfFluxFlag(), record0.getPsfFluxFlag()

