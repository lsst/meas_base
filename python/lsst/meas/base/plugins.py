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
"""
New plugins written purely in Python, based on meas_base plugin model
Plugins written purely in python.  Basis for the unit tests for this module,
and also a fall-back when no other Centroid algorithm is available.
"""
import numpy

from lsst.pex.exceptions import LsstCppException, LengthErrorException
from lsst.afw.table import tableLib
import lsst.afw.detection
import lsst.meas.algorithms

from .base import *
from .baseLib import *
from .sfm import *
from .forcedImage import *

WrappedSingleFramePlugin.generate(PsfFluxAlgorithm)
WrappedSingleFramePlugin.generate(SdssShapeAlgorithm)

class SingleFramePeakCentroidConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class PeakCentroid(SingleFramePlugin):
    """
    This is an easy to implement centroid algorithm, based on the peak pixel of the source
    footprint.  It is not the best centroid estimatation, but it is always available.
    """
    ConfigClass = SingleFramePeakCentroidConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        schema.addField("centroid.peak", "PointD", doc="measured centroid", units="pixels")
        schema.addField("centroid.peak.cov", "CovPointF", doc="covariance of measured centroid",
            units="pixels^2")

    def measure(self, exposure, source):
        peak = source.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        schema = source.getSchema()
        key = schema.find("centroid.peak").key
        source.set(key, result)

    def measureN(self, exposure, sources):
        raise NotImplementedError()

# plugin class must be registered to the singleton Registry class of the same type
# prior to any constructor requests (these are done in the __init__ of the measurement Task.
# will be executed during module initialization.

SingleFramePlugin.registry.register("centroid.peak", PeakCentroid)

# --- Forced Plugins ---

class ForcedPeakCentroidConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class ForcedPeakCentroid(ForcedPlugin):
    """
    The forced peak centroid is like the sfm peak centroid plugin, except that it must transform
    the peak coordinate from the original (reference) coordinate system to the coordinate system
    of the exposure being measured.
    """
    ConfigClass = ForcedPeakCentroidConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        field = tableLib.Field_PointD("centroid.peak", "measured peak centroid", "pixels")
        self.peakkey = schemaMapper.addOutputField(field)
        field = tableLib.Field_CovPointF("centroid.peak.cov", "measured peak centroid covariance", "pixels")
        self.covkey = schemaMapper.addOutputField(field)

    def measure(self, exposure, source, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()
        peak = refRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not referenceWcs == targetWcs:
            result = targetWcs.skyToPixel(referenceWcs.pixelToSky(result))
        source.set(self.peakkey, result)

    def measureN(self, exposure, measCat, refCat, refWcs):
        raise NotImplementedError()


ForcedPlugin.registry.register("centroid.peak", ForcedPeakCentroid)

class ForcedTransformedCentroidConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class ForcedTransformedCentroid(ForcedPlugin):

    ConfigClass = ForcedTransformedCentroidConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        field = tableLib.Field_PointD("transformed.peak", "measured peak centroid", "pixels")
        self.peakkey = schemaMapper.addOutputField(field)
        field = tableLib.Field_CovPointF("transformed.peak.cov", "measured peak centroid covariance",
            "pixels")
        self.covkey = schemaMapper.addOutputField(field)

    def measure(self, exposure, source, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()
        footprint = refRecord.getFootprint()
        peak = footprint.getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not targetWcs == referenceWcs:
            result = targetWcs.skyToPixel(referenceWcs.pixelToSky(result))
        source.set(self.peakkey, result)

    def measureN(self, exposure, measCat, refCat, refWcs):
        raise NotImplementedError()


# plugin class must be registered prior the singleton Registry class of the same type
# prior to any constructor requests (these are done in the __init__ of the measurement Task.
# will be executed during module initialization.
ForcedPlugin.registry.register("centroid.transformed", ForcedTransformedCentroid)
