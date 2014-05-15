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


def compareArrays(array1, array2, relDiff):
    if not array1.shape[0] == array2.shape[0] or not array1.shape[1] == array2.shape[1]:
         return False
    for i in range(array1.shape[0]):
        for j in range(array1.shape[1]):
            val1 = array1[i][j]
            val2 = array2[i][j]
            if numpy.isnan(val1) and numpy.isnan(val2):
                continue
            if numpy.isnan(val1) or numpy.isnan(val2):
                return False
            if abs(val1 == val2): continue
            if abs(2.0*(val1-val2)/(val1+val2)) < relDiff: continue
            return False
    return True

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
    print measCat.getShapeDefinition()
    print measCat.getPsfFluxDefinition()
    print measCat.getModelFluxDefinition()
    print measCat0.getCentroidDefinition()
    print measCat0.getShapeDefinition()
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

        value = record.getShape()
        value0 = record0.getShape()
        error = record.getShapeErr()
        error0 = record0.getShapeErr()
        flag0 = record0.getShapeFlag()
        flag = record.getShapeFlag() or record.get("base_SdssShape_flag_unweightedBad") or record.get("base_SdssShape_flag_unweighted") or record.get("base_SdssShape_flag_maxIter") or record.get("base_SdssShape_flag_shift")
        label = "Shape: "
        if not (value==value0):
            if not (numpy.isnan(value.getIxx()) and numpy.isnan(value0.getIxx())):
                print "Shape Values: ", record.getId(), value, value0
        if not compareArrays(error,error0, .001):
            print "Shape Errors: ", record.getId(), record.getCentroid()
        if not (flag == flag0):
            print "Shape Flags: ", record.getId(), flag,value,flag0,value0 
            print record.getShapeFlag(), record.get("base_SdssShape_flag_unweightedBad"), record.get("base_SdssShape_flag_unweighted"), record.get("base_SdssShape_flag_maxIter"), record.get("base_SdssShape_flag_shift")
            print record0.getShapeFlag(), record0.get("shape.sdss.flags.unweightedbad"), record0.get("shape.sdss.flags.unweighted"), record0.get("shape.sdss.flags.maxiter"), record0.get("shape.sdss.flags.shift")
