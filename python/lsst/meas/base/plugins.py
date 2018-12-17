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

"""Definition of measurement plugins.

This module defines and registers a series of pure-Python measurement plugins
which have trivial implementations. It also wraps measurement algorithms
defined in C++ to expose them to the measurement framework.
"""

import numpy as np

import lsst.pex.exceptions
import lsst.geom
import lsst.afw.detection
import lsst.afw.geom

from .pluginRegistry import register
from .pluginsBase import BasePlugin
from .baseMeasurement import BaseMeasurementPluginConfig
from .sfm import SingleFramePluginConfig, SingleFramePlugin
from .forcedMeasurement import ForcedPluginConfig, ForcedPlugin
from .wrappers import wrapSimpleAlgorithm, wrapTransform, GenericPlugin
from .transforms import SimpleCentroidTransform

from .apertureFlux import ApertureFluxControl, ApertureFluxTransform
from .transform import BaseTransform
from .blendedness import BlendednessAlgorithm, BlendednessControl
from .circularApertureFlux import CircularApertureFluxAlgorithm
from .gaussianFlux import GaussianFluxAlgorithm, GaussianFluxControl, GaussianFluxTransform
from .exceptions import MeasurementError
from .localBackground import LocalBackgroundControl, LocalBackgroundAlgorithm, LocalBackgroundTransform
from .naiveCentroid import NaiveCentroidAlgorithm, NaiveCentroidControl, NaiveCentroidTransform
from .peakLikelihoodFlux import PeakLikelihoodFluxAlgorithm, PeakLikelihoodFluxControl, \
    PeakLikelihoodFluxTransform
from .pixelFlags import PixelFlagsAlgorithm, PixelFlagsControl
from .psfFlux import PsfFluxAlgorithm, PsfFluxControl, PsfFluxTransform
from .scaledApertureFlux import ScaledApertureFluxAlgorithm, ScaledApertureFluxControl, \
    ScaledApertureFluxTransform
from .sdssCentroid import SdssCentroidAlgorithm, SdssCentroidControl, SdssCentroidTransform
from .sdssShape import SdssShapeAlgorithm, SdssShapeControl, SdssShapeTransform

__all__ = (
    "SingleFrameFPPositionConfig", "SingleFrameFPPositionPlugin",
    "SingleFrameJacobianConfig", "SingleFrameJacobianPlugin",
    "VarianceConfig", "SingleFrameVariancePlugin", "ForcedVariancePlugin",
    "InputCountConfig", "SingleFrameInputCountPlugin", "ForcedInputCountPlugin",
    "SingleFramePeakCentroidConfig", "SingleFramePeakCentroidPlugin",
    "SingleFrameSkyCoordConfig", "SingleFrameSkyCoordPlugin",
    "ForcedPeakCentroidConfig", "ForcedPeakCentroidPlugin",
    "ForcedTransformedCentroidConfig", "ForcedTransformedCentroidPlugin",
    "ForcedTransformedShapeConfig", "ForcedTransformedShapePlugin",
)


wrapSimpleAlgorithm(PsfFluxAlgorithm, Control=PsfFluxControl,
                    TransformClass=PsfFluxTransform, executionOrder=BasePlugin.FLUX_ORDER,
                    shouldApCorr=True, hasLogName=True)
wrapSimpleAlgorithm(PeakLikelihoodFluxAlgorithm, Control=PeakLikelihoodFluxControl,
                    TransformClass=PeakLikelihoodFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)
wrapSimpleAlgorithm(GaussianFluxAlgorithm, Control=GaussianFluxControl,
                    TransformClass=GaussianFluxTransform, executionOrder=BasePlugin.FLUX_ORDER,
                    shouldApCorr=True)
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

wrapSimpleAlgorithm(LocalBackgroundAlgorithm, Control=LocalBackgroundControl,
                    TransformClass=LocalBackgroundTransform, executionOrder=BasePlugin.FLUX_ORDER)

wrapTransform(PsfFluxTransform)
wrapTransform(PeakLikelihoodFluxTransform)
wrapTransform(GaussianFluxTransform)
wrapTransform(NaiveCentroidTransform)
wrapTransform(SdssCentroidTransform)
wrapTransform(SdssShapeTransform)
wrapTransform(ScaledApertureFluxTransform)
wrapTransform(ApertureFluxTransform)
wrapTransform(LocalBackgroundTransform)


class SingleFrameFPPositionConfig(SingleFramePluginConfig):
    """Configuration for the focal plane position measurment algorithm.
    """

    pass


@register("base_FPPosition")
class SingleFrameFPPositionPlugin(SingleFramePlugin):
    """Algorithm to calculate the position of a centroid on the focal plane.

    Parameters
    ----------
    config : `SingleFrameFPPositionConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schema : `lsst.afw.table.Schema`
        The schema for the measurement output catalog. New fields will be
        added to hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

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
            fp = lsst.geom.Point2D(np.nan, np.nan)
        else:
            center = measRecord.getCentroid()
            fp = det.transform(center, lsst.afw.cameraGeom.PIXELS, lsst.afw.cameraGeom.FOCAL_PLANE)
        measRecord.set(self.focalValue, fp)

    def fail(self, measRecord, error=None):
        measRecord.set(self.focalFlag, True)


class SingleFrameJacobianConfig(SingleFramePluginConfig):
    """Configuration for the Jacobian calculation plugin.
    """

    pixelScale = lsst.pex.config.Field(dtype=float, default=0.5, doc="Nominal pixel size (arcsec)")


@register("base_Jacobian")
class SingleFrameJacobianPlugin(SingleFramePlugin):
    """Compute the Jacobian and its ratio with a nominal pixel area.

    This enables one to compare relative, rather than absolute, pixel areas.

    Parameters
    ----------
    config : `SingleFrameJacobianConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schema : `lsst.afw.table.Schema`
        The schema for the measurement output catalog. New fields will be
        added to hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

    ConfigClass = SingleFrameJacobianConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.jacValue = schema.addField(name + '_value', type="D", doc="Jacobian correction")
        self.jacFlag = schema.addField(name + '_flag', type="Flag", doc="Set to 1 for any fatal failure")
        # Calculate one over the area of a nominal reference pixel, where area is in arcsec^2
        self.scale = pow(self.config.pixelScale, -2)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()
        # Compute the area of a pixel at a source record's centroid, and take
        # the ratio of that with the defined reference pixel area.
        result = np.abs(self.scale*exposure.getWcs().linearizePixelToSky(
            center,
            lsst.geom.arcseconds).getLinear().computeDeterminant())
        measRecord.set(self.jacValue, result)

    def fail(self, measRecord, error=None):
        measRecord.set(self.jacFlag, True)


class VarianceConfig(BaseMeasurementPluginConfig):
    """Configuration for the variance calculation plugin.
    """
    scale = lsst.pex.config.Field(dtype=float, default=5.0, optional=True,
                                  doc="Scale factor to apply to shape for aperture")
    mask = lsst.pex.config.ListField(doc="Mask planes to ignore", dtype=str,
                                     default=["DETECTED", "DETECTED_NEGATIVE", "BAD", "SAT"])


class VariancePlugin(GenericPlugin):
    """Compute the median variance corresponding to a footprint.

    The aim here is to measure the background variance, rather than that of
    the object itself. In order to achieve this, the variance is calculated
    over an area scaled up from the shape of the input footprint.

    Parameters
    ----------
    config : `VarianceConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schema : `lsst.afw.table.Schema`
        The schema for the measurement output catalog. New fields will be
        added to hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

    ConfigClass = VarianceConfig

    FAILURE_BAD_CENTROID = 1
    """Denotes failures due to bad centroiding (`int`).
    """

    FAILURE_EMPTY_FOOTPRINT = 2
    """Denotes failures due to a lack of usable pixels (`int`).
    """

    @classmethod
    def getExecutionOrder(cls):
        return BasePlugin.FLUX_ORDER

    def __init__(self, config, name, schema, metadata):
        GenericPlugin.__init__(self, config, name, schema, metadata)
        self.varValue = schema.addField(name + '_value', type="D", doc="Variance at object position")
        self.emptyFootprintFlag = schema.addField(name + '_flag_emptyFootprint', type="Flag",
                                                  doc="Set to True when the footprint has no usable pixels")

        # Alias the badCentroid flag to that which is defined for the target
        # of the centroid slot.  We do not simply rely on the alias because
        # that could be changed post-measurement.
        schema.getAliasMap().set(name + '_flag_badCentroid', schema.getAliasMap().apply("slot_Centroid_flag"))

    def measure(self, measRecord, exposure, center):
        # Create an aperture and grow it by scale value defined in config to
        # ensure there are enough pixels around the object to get decent
        # statistics
        if not np.all(np.isfinite(measRecord.getCentroid())):
            raise MeasurementError("Bad centroid and/or shape", self.FAILURE_BAD_CENTROID)
        aperture = lsst.afw.geom.Ellipse(measRecord.getShape(), measRecord.getCentroid())
        aperture.scale(self.config.scale)
        ellipse = lsst.afw.geom.SpanSet.fromShape(aperture)
        foot = lsst.afw.detection.Footprint(ellipse)
        foot.clipTo(exposure.getBBox(lsst.afw.image.PARENT))
        # Filter out any pixels which have mask bits set corresponding to the
        # planes to be excluded (defined in config.mask)
        maskedImage = exposure.getMaskedImage()
        pixels = lsst.afw.detection.makeHeavyFootprint(foot, maskedImage)
        maskBits = maskedImage.getMask().getPlaneBitMask(self.config.mask)
        logicalMask = np.logical_not(pixels.getMaskArray() & maskBits)
        # Compute the median variance value for each pixel not excluded by the
        # mask and write the record.  Numpy median is used here instead of
        # afw.math makeStatistics because of an issue with data types being
        # passed into the C++ layer (DM-2379).
        if np.any(logicalMask):
            medVar = np.median(pixels.getVarianceArray()[logicalMask])
            measRecord.set(self.varValue, medVar)
        else:
            raise MeasurementError("Footprint empty, or all pixels are masked, can't compute median",
                                   self.FAILURE_EMPTY_FOOTPRINT)

    def fail(self, measRecord, error=None):
        # Check that we have an error object and that it is of type
        # MeasurementError
        if isinstance(error, MeasurementError):
            assert error.getFlagBit() in (self.FAILURE_BAD_CENTROID, self.FAILURE_EMPTY_FOOTPRINT)
            # FAILURE_BAD_CENTROID handled by alias to centroid record.
            if error.getFlagBit() == self.FAILURE_EMPTY_FOOTPRINT:
                measRecord.set(self.emptyFootprintFlag, True)
        measRecord.set(self.varValue, np.nan)
        GenericPlugin.fail(self, measRecord, error)


SingleFrameVariancePlugin = VariancePlugin.makeSingleFramePlugin("base_Variance")
"""Single-frame version of `VariancePlugin`.
"""

ForcedVariancePlugin = VariancePlugin.makeForcedPlugin("base_Variance")
"""Forced version of `VariancePlugin`.
"""

class InputCountConfig(BaseMeasurementPluginConfig):
    """Configuration for the input image counting plugin.
    """
    pass


class InputCountPlugin(GenericPlugin):
    """Count the number of input images which contributed to a a source.

    Parameters
    ----------
    config : `InputCountConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schema : `lsst.afw.table.Schema`
        The schema for the measurement output catalog. New fields will be
        added to hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog

    Notes
    -----
    Information is derived from the image's `~lsst.afw.image.CoaddInputs`.
    Note these limitation:

    - This records the number of images which contributed to the pixel in the
      center of the source footprint, rather than to any or all pixels in the
      source.
    - Clipping in the coadd is not taken into account.
    """

    ConfigClass = InputCountConfig

    FAILURE_BAD_CENTROID = 1
    """Denotes failures due to bad centroiding (`int`).
    """

    FAILURE_NO_INPUTS = 2
    """Denotes failures due to the image not having coadd inputs.  (`int`)
    """

    @classmethod
    def getExecutionOrder(cls):
        return BasePlugin.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        GenericPlugin.__init__(self, config, name, schema, metadata)
        self.numberKey = schema.addField(name + '_value', type="I",
                                         doc="Number of images contributing at center, not including any"
                                             "clipping")
        self.noInputsFlag = schema.addField(name + '_flag_noInputs', type="Flag",
                                            doc="No coadd inputs available")
        # Alias the badCentroid flag to that which is defined for the target of the centroid slot.
        # We do not simply rely on the alias because that could be changed post-measurement.
        schema.getAliasMap().set(name + '_flag_badCentroid', schema.getAliasMap().apply("slot_Centroid_flag"))

    def measure(self, measRecord, exposure, center):
        if not exposure.getInfo().getCoaddInputs():
            raise MeasurementError("No coadd inputs defined.", self.FAILURE_NO_INPUTS)
        if not np.all(np.isfinite(center)):
            raise MeasurementError("Source has a bad centroid.", self.FAILURE_BAD_CENTROID)

        ccds = exposure.getInfo().getCoaddInputs().ccds
        measRecord.set(self.numberKey, len(ccds.subsetContaining(center, exposure.getWcs())))

    def fail(self, measRecord, error=None):
        if error is not None:
            assert error.getFlagBit() in (self.FAILURE_BAD_CENTROID, self.FAILURE_NO_INPUTS)
            # FAILURE_BAD_CENTROID handled by alias to centroid record.
            if error.getFlagBit() == self.FAILURE_NO_INPUTS:
                measRecord.set(self.noInputsFlag, True)
        GenericPlugin.fail(self, measRecord, error)


SingleFrameInputCountPlugin = InputCountPlugin.makeSingleFramePlugin("base_InputCount")
"""Single-frame version of `InputCoutPlugin`.
"""

ForcedInputCountPlugin = InputCountPlugin.makeForcedPlugin("base_InputCount")
"""Forced version of `InputCoutPlugin`.
"""


class SingleFramePeakCentroidConfig(SingleFramePluginConfig):
    """Configuration for the single frame peak centroiding algorithm.
    """
    pass


@register("base_PeakCentroid")
class SingleFramePeakCentroidPlugin(SingleFramePlugin):
    """Record the highest peak in a source footprint as its centroid.

    This is of course a relatively poor measure of the true centroid of the
    object; this algorithm is provided mostly for testing and debugging.

    Parameters
    ----------
    config : `SingleFramePeakCentroidConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schema : `lsst.afw.table.Schema`
        The schema for the measurement output catalog. New fields will be
        added to hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
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
    """Configuration for the sky coordinates algorithm.
    """
    pass


@register("base_SkyCoord")
class SingleFrameSkyCoordPlugin(SingleFramePlugin):
    """Record the sky position of an object based on its centroid slot and WCS.

    The position is record in the ``coord`` field, which is part of the
    `~lsst.afw.table.SourceCatalog` minimal schema.

    Parameters
    ----------
    config : `SingleFrameSkyCoordConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schema : `lsst.afw.table.Schema`
        The schema for the measurement output catalog. New fields will be
        added to hold measurements produced by this plugin.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog
    """

    ConfigClass = SingleFrameSkyCoordConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def measure(self, measRecord, exposure):
        # There should be a base class method for handling this exception. Put
        # this on a later ticket. Also, there should be a python Exception of
        # the appropriate type for this error
        if not exposure.hasWcs():
            raise Exception("Wcs not attached to exposure.  Required for " + self.name + " algorithm")
        measRecord.updateCoord(exposure.getWcs())

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception: this is
        # not ideal, but we don't have a place to put failures because we
        # don't allocate any fields.  Should consider fixing as part of
        # DM-1011
        pass


class ForcedPeakCentroidConfig(ForcedPluginConfig):
    """Configuration for the forced peak centroid algorithm.
    """
    pass


@register("base_PeakCentroid")
class ForcedPeakCentroidPlugin(ForcedPlugin):
    """Record the highest peak in a source footprint as its centroid.

    This is of course a relatively poor measure of the true centroid of the
    object; this algorithm is provided mostly for testing and debugging.

    This is similar to `SingleFramePeakCentroidPlugin`, except that transforms
    the peak coordinate from the original (reference) coordinate system to the
    coordinate system of the exposure being measured.

    Parameters
    ----------
    config : `ForcedPeakCentroidConfig`
        Plugin configuraion.
    name : `str`
        Plugin name.
    schemaMapper : `lsst.afw.table.SchemaMapper`
        A mapping from reference catalog fields to output
        catalog fields. Output fields are added to the output schema.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog.
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
        result = lsst.geom.Point2D(peak.getFx(), peak.getFy())
        result = targetWcs.skyToPixel(refWcs.pixelToSky(result))
        measRecord.set(self.keyX, result.getX())
        measRecord.set(self.keyY, result.getY())

    @staticmethod
    def getTransformClass():
        return SimpleCentroidTransform


class ForcedTransformedCentroidConfig(ForcedPluginConfig):
    """Configuration for the forced transformed centroid algorithm.
    """
    pass


@register("base_TransformedCentroid")
class ForcedTransformedCentroidPlugin(ForcedPlugin):
    """Record the transformation of the reference catalog centroid.

    The centroid recorded in the reference catalog is tranformed to the
    measurement coordinate system and stored.

    Parameters
    ----------
    config : `ForcedTransformedCentroidConfig`
        Plugin configuration
    name : `str`
        Plugin name
    schemaMapper : `lsst.afw.table.SchemaMapper`
        A mapping from reference catalog fields to output
        catalog fields. Output fields are added to the output schema.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog.

    Notes
    -----
    This is used as the slot centroid by default in forced measurement,
    allowing subsequent measurements to simply refer to the slot value just as
    they would in single-frame measurement.
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
    """Configuration for the forced transformed shape algorithm.
    """
    pass


@register("base_TransformedShape")
class ForcedTransformedShapePlugin(ForcedPlugin):
    """Record the transformation of the reference catalog shape.

    The shape recorded in the reference catalog is tranformed to the
    measurement coordinate system and stored.

    Parameters
    ----------
    config : `ForcedTransformedShapeConfig`
        Plugin configuration
    name : `str`
        Plugin name
    schemaMapper : `lsst.afw.table.SchemaMapper`
        A mapping from reference catalog fields to output
        catalog fields. Output fields are added to the output schema.
    metadata : `lsst.daf.base.PropertySet`
        Plugin metadata that will be attached to the output catalog.

    Notes
    -----
    This is used as the slot shape by default in forced measurement, allowing
    subsequent measurements to simply refer to the slot value just as they
    would in single-frame measurement.
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
            fullTransform = lsst.afw.geom.makeWcsPairTransform(refWcs, targetWcs)
            localTransform = lsst.afw.geom.linearizeTransform(fullTransform, refRecord.getCentroid())
            measRecord.set(self.shapeKey, refRecord.getShape().transform(localTransform.getLinear()))
        else:
            measRecord.set(self.shapeKey, refRecord.getShape())
        if self.flagKey is not None:
            measRecord.set(self.flagKey, refRecord.getShapeFlag())
