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
Definitions and registration of pure-Python plugins with trivial implementations,
and automatic plugin-from-algorithm calls for those implemented in C++.
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

# --- Wrapped C++ Plugins ---

WrappedSingleFramePlugin.generate(PsfFluxAlgorithm)
WrappedSingleFramePlugin.generate(SdssShapeAlgorithm)
WrappedSingleFramePlugin.generate(ApertureFluxAlgorithm)
WrappedSingleFramePlugin.generate(NaiveCentroidAlgorithm)
#WrappedSingleFramePlugin.generate(SdssCentroidAlgorithm)

WrappedForcedPlugin.generate(PsfFluxAlgorithm)

# --- Single-Frame Measurement Plugins ---

class SingleFramePeakCentroidConfig(SingleFramePluginConfig):

    def setDefaults(self):
        SingleFramePluginConfig.setDefaults(self)
        self.executionOrder = 0.0

class SingleFramePeakCentroidPlugin(SingleFramePlugin):
    """
    This is an easy-to-implement centroid algorithm, based on the peak pixel of the source
    footprint.  It is not the best centroid estimatation, but it is always available.
    """
    ConfigClass = SingleFramePeakCentroidConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.key = schema.addField(name, type="PointD", doc="peak centroid", units="pixels")

    def measure(self, measRecord, exposure):
        peak = measRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        measRecord.set(self.key, result)

# Plugin class must be registered to the singleton Registry class of the same type in order to
# be available for use with the corresponding measurement Task.
SingleFramePlugin.registry.register("centroid.peak", SingleFramePeakCentroidPlugin)


# --- Forced Plugins ---

class ForcedPeakCentroidConfig(ForcedPluginConfig):

    def setDefaults(self):
        ForcedPluginConfig.setDefaults(self)
        self.executionOrder = 0.0

class ForcedPeakCentroidPlugin(ForcedPlugin):
    """
    The forced peak centroid is like the SFM peak centroid plugin, except that it must transform
    the peak coordinate from the original (reference) coordinate system to the coordinate system
    of the exposure being measured.
    """
    ConfigClass = ForcedPeakCentroidConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.editOutputSchema()
        self.key = schema.addField(name, type="PointD", doc="peak centroid", units="pixels")

    def measure(self, measRecord, exposure, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()
        peak = refRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not referenceWcs == targetWcs:
            result = targetWcs.skyToPixel(referenceWcs.pixelToSky(result))
        measRecord.set(self.key, result)

ForcedPlugin.registry.register("centroid.peak", ForcedPeakCentroidPlugin)


class ForcedTransformedCentroidConfig(ForcedPluginConfig):

    def setDefaults(self):
        ForcedPluginConfig.setDefaults(self)
        self.executionOrder = 0.0

class ForcedTransformedCentroidPlugin(ForcedPlugin):
    """A centroid "algorithm" for forced measurement that simply transforms the centroid
    from the reference catalog to the measurement coordinate system.  This is used as
    the slot centroid by default in forced measurement, allowing subsequent measurements
    to simply refer to the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedTransformedCentroidConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.editOutputSchema()
        self.key = schema.addField(name, type="PointD", doc="transformed reference centroid", units="pixels")

    def measure(self, measRecord, exposure, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()
        # Note: may be better to transform refRecord.getCentroid() all the way once slots are working better
        result = targetWcs.skyToPixel(refRecord.getCoord())
        measRecord.set(self.key, result)

ForcedPlugin.registry.register("centroid.transformed", ForcedTransformedCentroidPlugin)

