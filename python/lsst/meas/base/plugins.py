#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
import lsst.afw.detection
import lsst.afw.geom

from .pluginRegistry import register
from .pluginsBase import BasePlugin
from .sfm import SingleFramePluginConfig, SingleFramePlugin
from .forcedMeasurement import ForcedPluginConfig, ForcedPlugin
from .wrappers import wrapSimpleAlgorithm, wrapTransform
from .transforms import SimpleCentroidTransform

from .apertureFlux import ApertureFluxControl, ApertureFluxTransform
from .transform import BaseTransform
from .blendedness import BlendednessAlgorithm, BlendednessControl
from .circularApertureFlux import CircularApertureFluxAlgorithm
from .gaussianCentroid import GaussianCentroidAlgorithm, GaussianCentroidControl, GaussianCentroidTransform
from .gaussianFlux import GaussianFluxAlgorithm, GaussianFluxControl, GaussianFluxTransform
from .exceptions import MeasurementError
from .naiveCentroid import NaiveCentroidAlgorithm, NaiveCentroidControl, NaiveCentroidTransform
from .peakLikelihoodFlux import PeakLikelihoodFluxAlgorithm, PeakLikelihoodFluxControl, PeakLikelihoodFluxTransform
from .pixelFlags import PixelFlagsAlgorithm, PixelFlagsControl
from .psfFlux import PsfFluxAlgorithm, PsfFluxControl, PsfFluxTransform
from .scaledApertureFlux import ScaledApertureFluxAlgorithm, ScaledApertureFluxControl, ScaledApertureFluxTransform
from .sdssCentroid import SdssCentroidAlgorithm, SdssCentroidControl, SdssCentroidTransform
from .sdssShape import SdssShapeAlgorithm, SdssShapeControl, SdssShapeTransform

__all__ = (
    "SingleFrameFPPositionConfig", "SingleFrameFPPositionPlugin",
    "SingleFrameJacobianConfig", "SingleFrameJacobianPlugin",
    "SingleFrameVarianceConfig", "SingleFrameVariancePlugin",
    "SingleFrameInputCountConfig", "SingleFrameInputCountPlugin",
    "SingleFramePeakCentroidConfig", "SingleFramePeakCentroidPlugin",
    "SingleFrameSkyCoordConfig", "SingleFrameSkyCoordPlugin",
    "ForcedPeakCentroidConfig", "ForcedPeakCentroidPlugin",
    "ForcedTransformedCentroidConfig", "ForcedTransformedCentroidPlugin",
    "ForcedTransformedShapeConfig", "ForcedTransformedShapePlugin",
)

# --- Wrapped C++ Plugins ---

wrapSimpleAlgorithm(PsfFluxAlgorithm, Control=PsfFluxControl,
                    TransformClass=PsfFluxTransform, executionOrder=BasePlugin.FLUX_ORDER,
                    shouldApCorr=True, hasLogName=True)
wrapSimpleAlgorithm(PeakLikelihoodFluxAlgorithm, Control=PeakLikelihoodFluxControl,
                    TransformClass=PeakLikelihoodFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)
wrapSimpleAlgorithm(GaussianFluxAlgorithm, Control=GaussianFluxControl,
                    TransformClass=GaussianFluxTransform, executionOrder=BasePlugin.FLUX_ORDER,
                    shouldApCorr=True)
wrapSimpleAlgorithm(GaussianCentroidAlgorithm, Control=GaussianCentroidControl,
                    TransformClass=GaussianCentroidTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapSimpleAlgorithm(NaiveCentroidAlgorithm, Control=NaiveCentroidControl,
                    TransformClass=NaiveCentroidTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapSimpleAlgorithm(SdssCentroidAlgorithm, Control=SdssCentroidControl,
                    TransformClass=SdssCentroidTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapSimpleAlgorithm(PixelFlagsAlgorithm, Control=PixelFlagsControl,
                    executionOrder=BasePlugin.FLUX_ORDER)
wrapSimpleAlgorithm(SdssShapeAlgorithm, Control=SdssShapeControl,
                    TransformClass=SdssShapeTransform, executionOrder=BasePlugin.SHAPE_ORDER)
wrapSimpleAlgorithm(ScaledApertureFluxAlgorithm, Control=ScaledApertureFluxControl,
                    TransformClass=ScaledApertureFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)

wrapSimpleAlgorithm(CircularApertureFluxAlgorithm, needsMetadata=True, Control=ApertureFluxControl,
                    TransformClass=ApertureFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)
wrapSimpleAlgorithm(BlendednessAlgorithm, Control=BlendednessControl,
                    TransformClass=BaseTransform, executionOrder=BasePlugin.SHAPE_ORDER)

wrapTransform(PsfFluxTransform)
wrapTransform(PeakLikelihoodFluxTransform)
wrapTransform(GaussianFluxTransform)
wrapTransform(GaussianCentroidTransform)
wrapTransform(NaiveCentroidTransform)
wrapTransform(SdssCentroidTransform)
wrapTransform(SdssShapeTransform)
wrapTransform(ScaledApertureFluxTransform)
wrapTransform(ApertureFluxTransform)

# --- Single-Frame Measurement Plugins ---


class SingleFrameFPPositionConfig(SingleFramePluginConfig):
    pass


@register("base_FPPosition")
class SingleFrameFPPositionPlugin(SingleFramePlugin):
    '''
    Algorithm to calculate the position of a centroid on the focal plane
    '''

    ConfigClass = SingleFrameFPPositionConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.focalValue = lsst.afw.table.Point2DKey.addFields(schema, name, "Position on the focal plane",
                                                              "mm")
        self.focalFlag = schema.addField(name + "_flag", type="Flag", doc="Set to True for any fatal failure")
        self.detectorFlag = schema.addField(name + "_missingDetector_flag", type="Flag",
                                            doc="Set to True if detector object is missing")

    def measure(self, measRecord, exposure):
        det = exposure.getDetector()
        if not det:
            measRecord.set(self.detectorFlag, True)
            fp = lsst.afw.geom.Point2D(numpy.nan, numpy.nan)
        else:
            center = measRecord.getCentroid()
            posInPix = det.makeCameraPoint(center, lsst.afw.cameraGeom.PIXELS)
            fp = det.transform(posInPix, lsst.afw.cameraGeom.FOCAL_PLANE).getPoint()
        measRecord.set(self.focalValue, fp)

    def fail(self, measRecord, error=None):
        measRecord.set(self.focalFlag, True)


class SingleFrameJacobianConfig(SingleFramePluginConfig):
    pixelScale = lsst.pex.config.Field(dtype=float, default=0.5, doc="Nominal pixel size (arcsec)")


@register("base_Jacobian")
class SingleFrameJacobianPlugin(SingleFramePlugin):
    '''
    Algorithm which computes the Jacobian about a source and computes its ratio with a nominal pixel area.
    This allows one to compare relative instead of absolute areas of pixels.
    '''

    ConfigClass = SingleFrameJacobianConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.jacValue = schema.addField(name + '_value', type="D", doc="Jacobian correction")
        self.jacFlag = schema.addField(name + '_flag', type="Flag", doc="Set to 1 for any fatal failure")
        # Calculate one over the area of a nominal reference pixel
        self.scale = pow(self.config.pixelScale, -2)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()
        # Compute the area of a pixel at the center of a source records centroid, and take the
        # ratio of that with the defined reference pixel area.
        result = numpy.abs(self.scale*exposure.getWcs().linearizePixelToSky(
            center,
            lsst.afw.geom.arcseconds).getLinear().computeDeterminant())
        measRecord.set(self.jacValue, result)

    def fail(self, measRecord, error=None):
        measRecord.set(self.jacFlag, True)


class SingleFrameVarianceConfig(SingleFramePluginConfig):
    scale = lsst.pex.config.Field(dtype=float, default=5.0, optional=True,
                                  doc="Scale factor to apply to shape for aperture")
    mask = lsst.pex.config.ListField(doc="Mask planes to ignore", dtype=str,
                                     default=["DETECTED", "DETECTED_NEGATIVE", "BAD", "SAT"])


@register("base_Variance")
class SingleFrameVariancePlugin(SingleFramePlugin):
    '''
    Calculate the median variance within a Footprint scaled from the object shape so
    the value is not terribly influenced by the object and instead represents the
    variance in the background near the object.
    '''
    ConfigClass = SingleFrameVarianceConfig
    FAILURE_BAD_CENTROID = 1
    FAILURE_EMPTY_FOOTPRINT = 2

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.varValue = schema.addField(name + '_value', type="D", doc="Variance at object position")
        self.varFlag = schema.addField(name + '_flag', type="Flag", doc="Set to True for any fatal failure")
        self.emptyFootprintFlag = schema.addField(name + '_flag_emptyFootprint', type="Flag",
                                                  doc="Set to True when the footprint has no usable pixels")

        # Alias the badCentroid flag to that which is defined for the target of the centroid slot.
        # We do not simply rely on the alias because that could be changed post-measurement.
        schema.getAliasMap().set(name + '_flag_badCentroid', schema.getAliasMap().apply("slot_Centroid_flag"))

    def measure(self, measRecord, exposure):
        if measRecord.getCentroidFlag():
            raise MeasurementError("Source record has a bad centroid.", self.FAILURE_BAD_CENTROID)
        # Create an aperture and grow it by scale value defined in config to ensure there are enough
        # pixels around the object to get decent statistics
        aperture = lsst.afw.geom.ellipses.Ellipse(measRecord.getShape(), measRecord.getCentroid())
        aperture.scale(self.config.scale)
        ellipse = lsst.afw.geom.SpanSet.fromShape(aperture)
        foot = lsst.afw.detection.Footprint(ellipse)
        foot.clipTo(exposure.getBBox(lsst.afw.image.PARENT))
        # Filter out any pixels which have mask bits set corresponding to the planes to be excluded
        # (defined in config.mask)
        maskedImage = exposure.getMaskedImage()
        pixels = lsst.afw.detection.makeHeavyFootprint(foot, maskedImage)
        maskBits = maskedImage.getMask().getPlaneBitMask(self.config.mask)
        logicalMask = numpy.logical_not(pixels.getMaskArray() & maskBits)
        # Compute the median variance value for each pixel not excluded by the mask and write the record.
        # Numpy median is used here instead of afw.math makeStatistics because of an issue with data types
        # being passed into the C++ layer (DM-2379).
        if numpy.any(logicalMask):
            medVar = numpy.median(pixels.getVarianceArray()[logicalMask])
            measRecord.set(self.varValue, medVar)
        else:
            raise MeasurementError("Footprint empty, or all pixels are masked, can't compute median",
                                      self.FAILURE_EMPTY_FOOTPRINT)

    def fail(self, measRecord, error=None):
        # Check that we have a error object and that it is of type MeasurementError
        if isinstance(error, MeasurementError):
            assert error.getFlagBit() in (self.FAILURE_BAD_CENTROID, self.FAILURE_EMPTY_FOOTPRINT)
            # FAILURE_BAD_CENTROID handled by alias to centroid record.
            if error.getFlagBit() == self.FAILURE_EMPTY_FOOTPRINT:
                measRecord.set(self.emptyFootprintFlag, True)
        measRecord.set(self.varValue, numpy.nan)
        measRecord.set(self.varFlag, True)


class SingleFrameInputCountConfig(SingleFramePluginConfig):
    pass


@register("base_InputCount")
class SingleFrameInputCountPlugin(SingleFramePlugin):
    """
    Plugin to count how many input images contributed to each source. This information
    is in the exposure's coaddInputs. Some limitations:
    * This is only for the pixel containing the center, not for all the pixels in the
      Footprint
    * This does not account for any clipping in the coadd
    """

    ConfigClass = SingleFrameInputCountConfig
    FAILURE_BAD_CENTROID = 1
    FAILURE_NO_INPUTS = 2

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.numberKey = schema.addField(name + '_value', type="I",
                                         doc="Number of images contributing at center, not including any"
                                             "clipping")
        self.numberFlag = schema.addField(name + '_flag', type="Flag", doc="Set to True for fatal failure")
        self.noInputsFlag = schema.addField(name + '_flag_noInputs', type="Flag",
                                            doc="No coadd inputs available")
        # Alias the badCentroid flag to that which is defined for the target of the centroid slot.
        # We do not simply rely on the alias because that could be changed post-measurement.
        schema.getAliasMap().set(name + '_flag_badCentroid', schema.getAliasMap().apply("slot_Centroid_flag"))

    def measure(self, measRecord, exposure):
        if not exposure.getInfo().getCoaddInputs():
            raise MeasurementError("No coadd inputs defined.", self.FAILURE_NO_INPUTS)
        if measRecord.getCentroidFlag():
            raise MeasurementError("Source has a bad centroid.", self.FAILURE_BAD_CENTROID)

        center = measRecord.getCentroid()
        ccds = exposure.getInfo().getCoaddInputs().ccds
        measRecord.set(self.numberKey, len(ccds.subsetContaining(center, exposure.getWcs())))

    def fail(self, measRecord, error=None):
        measRecord.set(self.numberFlag, True)
        if error is not None:
            assert error.getFlagBit() in (self.FAILURE_BAD_CENTROID, self.FAILURE_NO_INPUTS)
            # FAILURE_BAD_CENTROID handled by alias to centroid record.
            if error.getFlagBit() == self.FAILURE_NO_INPUTS:
                measRecord.set(self.noInputsFlag, True)


class SingleFramePeakCentroidConfig(SingleFramePluginConfig):
    pass


@register("base_PeakCentroid")
class SingleFramePeakCentroidPlugin(SingleFramePlugin):
    """
    A centroid algorithm that simply uses the first (i.e. highest) Peak in the Source's
    Footprint as the centroid.  This is of course a relatively poor measure of the true
    centroid of the object; this algorithm is provided mostly for testing and debugging.
    """

    ConfigClass = SingleFramePeakCentroidConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixel")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixel")
        self.flag = schema.addField(name + "_flag", type="Flag", doc="Centroiding failed")

    def measure(self, measRecord, exposure):
        peak = measRecord.getFootprint().getPeaks()[0]
        measRecord.set(self.keyX, peak.getFx())
        measRecord.set(self.keyY, peak.getFy())

    def fail(self, measRecord, error=None):
        measRecord.set(self.flag, True)

    @staticmethod
    def getTransformClass():
        return SimpleCentroidTransform


class SingleFrameSkyCoordConfig(SingleFramePluginConfig):
    pass


@register("base_SkyCoord")
class SingleFrameSkyCoordPlugin(SingleFramePlugin):
    """
    A measurement plugin that sets the "coord" field (part of the Source minimal schema)
    using the slot centroid and the Wcs attached to the Exposure.
    """

    ConfigClass = SingleFrameSkyCoordConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

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


# --- Forced Plugins ---

class ForcedPeakCentroidConfig(ForcedPluginConfig):
    pass


@register("base_PeakCentroid")
class ForcedPeakCentroidPlugin(ForcedPlugin):
    """
    The forced peak centroid is like the SFM peak centroid plugin, except that it must transform
    the peak coordinate from the original (reference) coordinate system to the coordinate system
    of the exposure being measured.
    """

    ConfigClass = ForcedPeakCentroidConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixel")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixel")

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        peak = refRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not refWcs == targetWcs:
            result = targetWcs.skyToPixel(refWcs.pixelToSky(result))
        measRecord.set(self.keyX, result.getX())
        measRecord.set(self.keyY, result.getY())

    @staticmethod
    def getTransformClass():
        return SimpleCentroidTransform


class ForcedTransformedCentroidConfig(ForcedPluginConfig):
    pass


@register("base_TransformedCentroid")
class ForcedTransformedCentroidPlugin(ForcedPlugin):
    """A centroid pseudo-algorithm for forced measurement that simply transforms the centroid
    from the reference catalog to the measurement coordinate system.  This is used as
    the slot centroid by default in forced measurement, allowing subsequent measurements
    to simply refer to the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedTransformedCentroidConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate x and y fields, join these into a single FunctorKey for ease-of-use.
        xKey = schema.addField(name + "_x", type="D", doc="transformed reference centroid column",
                               units="pixel")
        yKey = schema.addField(name + "_y", type="D", doc="transformed reference centroid row",
                               units="pixel")
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
        else:
            measRecord.set(self.centroidKey, refRecord.getCentroid())
        if self.flagKey is not None:
            measRecord.set(self.flagKey, refRecord.getCentroidFlag())


class ForcedTransformedShapeConfig(ForcedPluginConfig):
    pass


@register("base_TransformedShape")
class ForcedTransformedShapePlugin(ForcedPlugin):
    """A shape pseudo-algorithm for forced measurement that simply transforms the shape
    from the reference catalog to the measurement coordinate system.  This is used as
    the slot shape by default in forced measurement, allowing subsequent measurements
    to simply refer to the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedTransformedShapeConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schemaMapper, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate xx, yy, xy fields, join these into a single FunctorKey for ease-of-use.
        xxKey = schema.addField(name + "_xx", type="D", doc="transformed reference shape x^2 moment",
                                units="pixel^2")
        yyKey = schema.addField(name + "_yy", type="D", doc="transformed reference shape y^2 moment",
                                units="pixel^2")
        xyKey = schema.addField(name + "_xy", type="D", doc="transformed reference shape xy moment",
                                units="pixel^2")
        self.shapeKey = lsst.afw.table.QuadrupoleKey(xxKey, yyKey, xyKey)
        # Because we're taking the reference position as given, we don't bother transforming its
        # uncertainty and reporting that here, so there are no sigma or cov fields.  We do propagate
        # the flag field, if it exists.
        if "slot_Shape_flag" in schemaMapper.getInputSchema():
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
        else:
            measRecord.set(self.shapeKey, refRecord.getShape())
        if self.flagKey is not None:
            measRecord.set(self.flagKey, refRecord.getShapeFlag())
