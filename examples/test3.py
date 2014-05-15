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
import sys 
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

#  Read a catalog for the old (meas_algorithm) and new (meas_base) algorithms
#  These are the result of a complete run of processCcd with default configurations of measurement.py
#  and sfm.py.  "0" means the old meas_algorithm algorithms.  Only the measurement task is different.
#  The pipeline preparatory to measurement should be all the same for both catalogs.   

    errorLimit = 1
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
        value = record.get("classification_extendedness")
        value0= record0.get("classification.extendedness")
        label = "Classification: "
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlags: "
        value = record.get("base_PixelFlags_flag_edge")
        value0 = record0.get("flags.pixel.edge")
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlagsEdge: "
        value = record.get("base_PixelFlags_flag_interpolated")
        value0 = record0.get("flags.pixel.interpolated.any")
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlagPixelInterpolated: "
        value = record.get("base_PixelFlags_flag_interpolatedCenter")
        value0 = record0.get("flags.pixel.interpolated.center")
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlagInterpolatedCenter: "
        value = record.get("base_PixelFlags_flag_saturated")
        value0 = record0.get("flags.pixel.saturated.any")
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlagSaturated: "
        value = record.get("base_PixelFlags_flag_saturatedCenter")
        value0 = record0.get("flags.pixel.saturated.center")
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlagsCr: "
        value = record.get("base_PixelFlags_flag_cr")
        value0 = record0.get("flags.pixel.cr.any")
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        label = "PixelFlagsCrCenter: "
        value = record.get("base_PixelFlags_flag_crCenter")
        value0 = record0.get("flags.pixel.cr.center")
        label = "PixelFlagBad: "
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0
        value = record.get("base_PixelFlags_flag_bad")
        value0 = record0.get("flags.pixel.bad")
        label = "PixelFlags: "
        if not (value == value0) and  not(numpy.isnan(value)):
            print label, "Values differ: ", record.getCentroid(), record.getId(), value, value0

