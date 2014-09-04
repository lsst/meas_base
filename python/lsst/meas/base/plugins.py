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

import lsst.pex.exceptions
from lsst.afw.table import tableLib
import lsst.afw.detection

from .base import *
from .baseLib import *
from .sfm import *
from .forcedMeasurement import *

# --- Wrapped C++ Plugins ---

WrappedSingleFramePlugin.generate(PsfFluxAlgorithm)
WrappedSingleFramePlugin.generate(SdssShapeAlgorithm, executionOrder=1.0)
WrappedSingleFramePlugin.generate(SdssCentroidAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(SincFluxAlgorithm)
WrappedSingleFramePlugin.generate(PixelFlagsAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(NaiveFluxAlgorithm)
WrappedSingleFramePlugin.generate(GaussianCentroidAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(GaussianFluxAlgorithm)
WrappedSingleFramePlugin.generate(NaiveCentroidAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(PeakLikelihoodFluxAlgorithm)
WrappedForcedPlugin.generate(PsfFluxAlgorithm)

# --- Single-Frame Measurement Plugins ---

class SingleFramePeakCentroidConfig(SingleFramePluginConfig):

    def setDefaults(self):
        SingleFramePluginConfig.setDefaults(self)
        self.executionOrder = 0.0

class SingleFramePeakCentroidPlugin(SingleFramePlugin):
    """
    A centroid algorithm that simply uses the first (i.e. highest) Peak in the Source's
    Footprint as the centroid.  This is of course a relatively poor measure of the true
    centroid of the object; this algorithm is provided mostly for testing and debugging.
    """

    ConfigClass = SingleFramePeakCentroidConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixels")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixels")

    def measure(self, measRecord, exposure):
        peak = measRecord.getFootprint().getPeaks()[0]
        measRecord.set(self.keyX, peak.getFx())
        measRecord.set(self.keyY, peak.getFy())

SingleFramePlugin.registry.register("base_PeakCentroid", SingleFramePeakCentroidPlugin)


class SingleFrameSkyCoordConfig(SingleFramePluginConfig):

    usePeak = lsst.pex.config.Field(dtype=bool, default=False, optional=True,
                                  doc="use footprint peak instead of centroid slot ")
    def setDefaults(self):
        SingleFramePluginConfig.setDefaults(self)
        self.executionOrder = 5.0

class SingleFrameSkyCoordPlugin(SingleFramePlugin):
    """
    A measurement plugin that sets the "coord" field (part of the Source minimal schema)
    using the slot centroid and the Wcs attached to the Exposure.
    """
    ConfigClass = SingleFrameSkyCoordConfig

    def measure(self, measRecord, exposure):
        # there should be a base class method for handling this exception. Put this on a later ticket
        # Also, there should be a python Exception of the appropriate type for this error
        if not exposure.hasWcs():
            raise Exception("Wcs not attached to exposure.  Required for " + self.name + " algorithm")
        # this is a temporary hack around, since centroid don't work yet
        if self.config.usePeak:
            peak = measRecord.getFootprint().getPeaks()[0]
            coord = exposure.getWcs().pixelToSky(peak.getFx(), peak.getFy())
            measRecord.setCoord(coord)
        else:
            measRecord.updateCoord(exposure.getWcs())

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception: this is not ideal,
        # but we don't have a place to put failures because we don't allocate any fields.
        # Should consider fixing as part of DM-1011
        pass

SingleFramePlugin.registry.register("base_SkyCoord", SingleFrameSkyCoordPlugin)

class SingleFrameClassificationConfig(SingleFramePluginConfig):

    fluxRatio = lsst.pex.config.Field(dtype=float, default=.925, optional=True,
                                  doc="critical ratio of model to psf flux")
    modelErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                  doc="correction factor for modelFlux error")
    psfErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                  doc="correction factor for psfFlux error")
    def setDefaults(self):
        SingleFramePluginConfig.setDefaults(self)
        self.executionOrder = 5.0

class SingleFrameClassificationPlugin(SingleFramePlugin):
    """
    A binary measure of the extendedness of a source, based a simple cut on the ratio of the
    PSF flux to the model flux.

    Because the fluxes on which this algorithm is based are slot measurements, they can be provided
    by different algorithms, and the "fluxRatio" threshold used by this algorithm should generally
    be set differently for different algorithms.  To do this, plot the difference between the PSF
    magnitude and the model magnitude vs. the PSF magnitude, and look for where the cloud of galaxies
    begins.
    """
    ConfigClass = SingleFrameClassificationConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for extended sources, 0 for point sources.")

    def measure(self, measRecord, exposure):
        modelFlux = measRecord.getModelFlux()
        modelFluxErr = measRecord.getModelFluxErr()
        psfFlux = measRecord.getPsfFlux()
        psfFluxErr = measRecord.getPsfFluxErr()
        flux1 = self.config.fluxRatio*modelFlux + self.config.modelErrFactor*modelFluxErr
        flux2 = psfFlux + self.config.psfErrFactor*psfFluxErr
        if flux1 < flux2:
            measRecord.set(self.keyProbability, 0.0)
        else:
            measRecord.set(self.keyProbability, 1.0);

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception.  We should be setting a flag
        # instead.
        pass

SingleFramePlugin.registry.register("base_ClassificationExtendedness", SingleFrameClassificationPlugin)

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
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixels")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixels")

    def measure(self, measRecord, exposure, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()
        peak = refRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not referenceWcs == targetWcs:
            result = targetWcs.skyToPixel(referenceWcs.pixelToSky(result))
        measRecord.set(self.keyX, result.getX())
        measRecord.set(self.keyY, result.getY())

ForcedPlugin.registry.register("base_PeakCentroid", ForcedPeakCentroidPlugin)


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

ForcedPlugin.registry.register("base_TransformedCentroid", ForcedTransformedCentroidPlugin)

