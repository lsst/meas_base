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
from .wrappers import *

# --- Wrapped C++ Plugins ---

WrappedSingleFramePlugin.generate(SdssShapeAlgorithm, executionOrder=1.0)
WrappedSingleFramePlugin.generate(SdssCentroidAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(SincFluxAlgorithm)
WrappedSingleFramePlugin.generate(PixelFlagsAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(NaiveFluxAlgorithm)
WrappedSingleFramePlugin.generate(GaussianCentroidAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(GaussianFluxAlgorithm)
WrappedSingleFramePlugin.generate(NaiveCentroidAlgorithm, executionOrder=0.0)
WrappedSingleFramePlugin.generate(PeakLikelihoodFluxAlgorithm)
WrappedForcedPlugin.generate(SdssShapeAlgorithm, executionOrder=1.0)
WrappedForcedPlugin.generate(SdssCentroidAlgorithm, executionOrder=0.0)
WrappedForcedPlugin.generate(SincFluxAlgorithm)
WrappedForcedPlugin.generate(PixelFlagsAlgorithm, executionOrder=0.0)
WrappedForcedPlugin.generate(NaiveFluxAlgorithm)
WrappedForcedPlugin.generate(GaussianCentroidAlgorithm, executionOrder=0.0)
WrappedForcedPlugin.generate(GaussianFluxAlgorithm)
WrappedForcedPlugin.generate(NaiveCentroidAlgorithm, executionOrder=0.0)
WrappedForcedPlugin.generate(PeakLikelihoodFluxAlgorithm)

wrapSimpleAlgorithm(PsfFluxAlgorithm, Control=PsfFluxControl)

# --- Aperture Flux Measurement Plugins ---

# The aperture flux algorithms are wrapped differently (and more verbosely) than the rest, because
# the C++ wrapping mechanism we originally came up with wasn't sufficiently general to handle
# algorithms with a dynamically-determined set of flag fields.  That will be fixed in the future
# in DM-1130.

class CircularApertureFluxSingleFramePlugin(SingleFramePlugin):
    """!
    Measure a sequence of circular aperture fluxes.

    See the C++ CircularApertureFluxAlgorithm class for more information.
    """

    ConfigClass = lsst.pex.config.makeConfigClass(
        ApertureFluxControl,
        base=SingleFramePlugin.ConfigClass
    )

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        for radius in self.config.radii:
            metadata.add("base_CircularApertureFlux_radii", radius)
        self.algorithm = CircularApertureFluxAlgorithm(config.makeControl(), name, schema)

    def measure(self, measRecord, exposure):
        self.algorithm.measure(measRecord, exposure)

SingleFramePlugin.registry.register("base_CircularApertureFlux", CircularApertureFluxSingleFramePlugin)

class CircularApertureFluxForcedPlugin(ForcedPlugin):
    """!
    Measure a sequence of circular aperture fluxes.

    See the C++ CircularApertureFluxAlgorithm class for more information.
    """

    ConfigClass = lsst.pex.config.makeConfigClass(
        ApertureFluxControl,
        base=ForcedPlugin.ConfigClass
    )

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        for radius in self.config.radii:
            metadata.add("base_CircularApertureFlux_radii", radius)
        schema = schemaMapper.editOutputSchema()
        self.algorithm = CircularApertureFluxAlgorithm(config.makeControl(), name, schema)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        self.algorithm.measure(measRecord, exposure)

ForcedPlugin.registry.register("base_CircularApertureFlux", CircularApertureFluxForcedPlugin)

# --- Single-Frame Measurement Plugins ---

class SingleFramePeakCentroidConfig(SingleFramePluginConfig):

    def setDefaults(self):
        SingleFramePluginConfig.setDefaults(self)
        self.executionOrder = 0.0

@register("base_PeakCentroid")
class SingleFramePeakCentroidPlugin(SingleFramePlugin):
    """
    A centroid algorithm that simply uses the first (i.e. highest) Peak in the Source's
    Footprint as the centroid.  This is of course a relatively poor measure of the true
    centroid of the object; this algorithm is provided mostly for testing and debugging.
    """

    ConfigClass = SingleFramePeakCentroidConfig

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixels")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixels")

    def measure(self, measRecord, exposure):
        peak = measRecord.getFootprint().getPeaks()[0]
        measRecord.set(self.keyX, peak.getFx())
        measRecord.set(self.keyY, peak.getFy())


class SingleFrameSkyCoordConfig(SingleFramePluginConfig):

    def setDefaults(self):
        SingleFramePluginConfig.setDefaults(self)
        self.executionOrder = 5.0

@register("base_SkyCoord")
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
        measRecord.updateCoord(exposure.getWcs())

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception: this is not ideal,
        # but we don't have a place to put failures because we don't allocate any fields.
        # Should consider fixing as part of DM-1011
        pass

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

@register("base_ClassificationExtendedness")
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

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
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


# --- Forced Plugins ---

class ForcedPeakCentroidConfig(ForcedPluginConfig):

    def setDefaults(self):
        ForcedPluginConfig.setDefaults(self)
        self.executionOrder = 0.0

@register("base_PeakCentroid")
class ForcedPeakCentroidPlugin(ForcedPlugin):
    """
    The forced peak centroid is like the SFM peak centroid plugin, except that it must transform
    the peak coordinate from the original (reference) coordinate system to the coordinate system
    of the exposure being measured.
    """
    ConfigClass = ForcedPeakCentroidConfig

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixels")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixels")

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        peak = refRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not refWcs == targetWcs:
            result = targetWcs.skyToPixel(refWcs.pixelToSky(result))
        measRecord.set(self.keyX, result.getX())
        measRecord.set(self.keyY, result.getY())


class ForcedTransformedCentroidConfig(ForcedPluginConfig):

    def setDefaults(self):
        ForcedPluginConfig.setDefaults(self)
        self.executionOrder = 0.0

@register("base_TransformedCentroid")
class ForcedTransformedCentroidPlugin(ForcedPlugin):
    """A centroid pseudo-algorithm for forced measurement that simply transforms the centroid
    from the reference catalog to the measurement coordinate system.  This is used as
    the slot centroid by default in forced measurement, allowing subsequent measurements
    to simply refer to the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedTransformedCentroidConfig

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate x and y fields, join these into a single FunctorKey for ease-of-use.
        xKey = schema.addField(name + "_x", type="D", doc="transformed reference centroid column",
                               units="pixels")
        yKey = schema.addField(name + "_y", type="D", doc="transformed reference centroid row",
                               units="pixels")
        self.centroidKey = lsst.afw.table.Point2DKey(xKey, yKey)
        # Because we're taking the reference position as given, we don't bother transforming its
        # uncertainty and reporting that here, so there are no sigma or cov fields.  We do propagate
        # the flag field, if it exists.
        if "slot_Centroid_flag" in schemaMapper.getInputSchema():
            self.flagKey = schema.addField(name + "_flag", type="Flag",
                                           doc="whether the reference centroid is marked as bad")
        else:
            self.flagKey = None

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        if not refWcs == targetWcs:
            targetPos = targetWcs.skyToPixel(refWcs.pixelToSky(refRecord.getCentroid()))
            measRecord.set(self.centroidKey, targetPos)
        if self.flagKey is not None:
            measRecord.set(self.flagKey, refRecord.getCentroidFlag())


class ForcedTransformedShapeConfig(ForcedPluginConfig):

    def setDefaults(self):
        ForcedPluginConfig.setDefaults(self)
        self.executionOrder = 1.0

@register("base_TransformedShape")
class ForcedTransformedShapePlugin(ForcedPlugin):
    """A shape pseudo-algorithm for forced measurement that simply transforms the shape
    from the reference catalog to the measurement coordinate system.  This is used as
    the slot shape by default in forced measurement, allowing subsequent measurements
    to simply refer to the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedTransformedShapeConfig

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate xx, yy, xy fields, join these into a single FunctorKey for ease-of-use.
        xxKey = schema.addField(name + "_xx", type="D", doc="transformed reference shape x^2 moment",
                                units="pixels^2")
        yyKey = schema.addField(name + "_yy", type="D", doc="transformed reference shape y^2 moment",
                                units="pixels^2")
        xyKey = schema.addField(name + "_xy", type="D", doc="transformed reference shape xy moment",
                                units="pixels^2")
        self.shapeKey = lsst.afw.table.QuadrupoleKey(xxKey, yyKey, xyKey)
        # Because we're taking the reference position as given, we don't bother transforming its
        # uncertainty and reporting that here, so there are no sigma or cov fields.  We do propagate
        # the flag field, if it exists.
        if "slot_shape_flag" in schemaMapper.getInputSchema():
            self.flagKey = schema.addField(name + "_flag", type="Flag",
                                           doc="whether the reference shape is marked as bad")
        else:
            self.flagKey = None

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        if not refWcs == targetWcs:
            fullTransform = lsst.afw.image.XYTransformFromWcsPair(targetWcs, refWcs)
            localTransform = fullTransform.linearizeForwardTransform(refRecord.getCentroid())
            measRecord.set(self.shapeKey, refRecord.getShape().transform(localTransform.getLinear()))
        if self.flagKey is not None:
            measRecord.set(self.flagKey, refRecord.getShapeFlag())
