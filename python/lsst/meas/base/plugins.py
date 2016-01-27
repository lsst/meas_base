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
import lsst.afw.math as afwMath

from .pluginRegistry import register
from . import baseLib as bl
from .baseMeasurement import BasePlugin
from .sfm import SingleFramePluginConfig, SingleFramePlugin
from .forcedMeasurement import ForcedPluginConfig, ForcedPlugin
from .wrappers import wrapSimpleAlgorithm
from .transforms import SimpleCentroidTransform

__all__ = (
    "SingleFrameFPPositionConfig", "SingleFrameFPPositionPlugin",
    "SingleFrameJacobianConfig", "SingleFrameJacobianPlugin",
    "SingleFrameVarianceConfig", "SingleFrameVariancePlugin",
    "SingleFrameInputCountConfig", "SingleFrameInputCountPlugin",
    "SingleFramePeakCentroidConfig", "SingleFramePeakCentroidPlugin", 
    "SingleFrameSkyCoordConfig", "SingleFrameSkyCoordPlugin",
    "SingleFrameClassificationConfig", "SingleFrameClassificationPlugin",
    "ForcedPeakCentroidConfig", "ForcedPeakCentroidPlugin",
    "ForcedTransformedCentroidConfig", "ForcedTransformedCentroidPlugin",
    "ForcedTransformedShapeConfig", "ForcedTransformedShapePlugin",
)

# --- Wrapped C++ Plugins ---

wrapSimpleAlgorithm(bl.PsfFluxAlgorithm, Control=bl.PsfFluxControl,
                TransformClass=bl.PsfFluxTransform, executionOrder=BasePlugin.FLUX_ORDER, shouldApCorr=True)
wrapSimpleAlgorithm(bl.PeakLikelihoodFluxAlgorithm, Control=bl.PeakLikelihoodFluxControl,
                TransformClass=bl.PeakLikelihoodFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)
wrapSimpleAlgorithm(bl.GaussianFluxAlgorithm, Control=bl.GaussianFluxControl,
                TransformClass=bl.GaussianFluxTransform, executionOrder=BasePlugin.FLUX_ORDER, shouldApCorr=True)
wrapSimpleAlgorithm(bl.GaussianCentroidAlgorithm, Control=bl.GaussianCentroidControl,
                TransformClass=bl.GaussianCentroidTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapSimpleAlgorithm(bl.NaiveCentroidAlgorithm, Control=bl.NaiveCentroidControl,
                TransformClass=bl.NaiveCentroidTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapSimpleAlgorithm(bl.SdssCentroidAlgorithm, Control=bl.SdssCentroidControl,
                TransformClass=bl.SdssCentroidTransform, executionOrder=BasePlugin.CENTROID_ORDER)
wrapSimpleAlgorithm(bl.PixelFlagsAlgorithm, Control=bl.PixelFlagsControl, executionOrder=BasePlugin.FLUX_ORDER)
wrapSimpleAlgorithm(bl.SdssShapeAlgorithm, Control=bl.SdssShapeControl,
                TransformClass=bl.SdssShapeTransform, executionOrder=BasePlugin.SHAPE_ORDER)
wrapSimpleAlgorithm(bl.ScaledApertureFluxAlgorithm, Control=bl.ScaledApertureFluxControl,
                TransformClass=bl.ScaledApertureFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)

wrapSimpleAlgorithm(bl.CircularApertureFluxAlgorithm, needsMetadata=True, Control=bl.ApertureFluxControl,
                    TransformClass=bl.ApertureFluxTransform, executionOrder=BasePlugin.FLUX_ORDER)

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
        result = numpy.abs(self.scale*exposure.getWcs().linearizePixelToSky(center,
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

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.varValue = schema.addField(name + '_value', type="D", doc="Variance at object position")
        self.varFlag = schema.addField(name + '_flag', type="Flag", doc="Set to True for any fatal failure")

        # Alias the badCentroid flag to that which is defined for the target of the centroid slot.
        # We do not simply rely on the alias because that could be changed post-measurement.
        schema.getAliasMap().set(name + '_flag_badCentroid', schema.getAliasMap().apply("slot_Centroid_flag"))

    def measure(self, measRecord, exposure):
        if measRecord.getCentroidFlag():
            raise bl.MeasurementError("Source record has a bad centroid.", self.FAILURE_BAD_CENTROID)
        center = measRecord.getCentroid()
        # Create an aperture and grow it by scale value defined in config to ensure there are enough
        # pixels around the object to get decent statistics
        aperture = lsst.afw.geom.ellipses.Ellipse(measRecord.getShape(), measRecord.getCentroid())
        aperture.scale(self.config.scale)
        foot = lsst.afw.detection.Footprint(aperture)
        foot.clipTo(exposure.getBBox(lsst.afw.image.PARENT))
        # Filter out any pixels which have mask bits set corresponding to the planes to be excluded
        # (defined in config.mask)
        maskedImage = exposure.getMaskedImage()
        maskBits = maskedImage.getMask().getPlaneBitMask(self.config.mask)
        logicalMask = numpy.logical_not(maskedImage.getMask().getArray() & maskBits)
        # Compute the median variance value for each pixel not excluded by the mask and write the record.
        # Numpy median is used here instead of afw.math makeStatistics because of an issue with data types
        # being passed into the C++ layer (DM-2379).
        medVar = numpy.median(maskedImage.getVariance().getArray()[logicalMask])
        measRecord.set(self.varValue, medVar)

    def fail(self, measRecord, error=None):
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

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.numberKey = schema.addField(name + '_value', type="I",
                                         doc="Number of images contributing at center, not including any"
                                             "clipping")
        self.numberFlag = schema.addField(name + '_flag', type="Flag", doc="Set to True for fatal failure")

    def measure(self, measRecord, exposure):
        if not exposure.getInfo().getCoaddInputs():
            raise lsst.pex.exceptions.RuntimeError("No coadd inputs defined")

        center = measRecord.getCentroid()
        # Promote bounding box of the Footprint to type D to ensure
        # the ability to compare the Footprint and center (they may be of mixed types I and D)
        fpBbox = lsst.afw.geom.Box2D(measRecord.getFootprint().getBBox())

        # Check to ensure that the center exists and that it is contained within the Footprint
        # to catch bad centroiding
        if not center:
            raise Exception("The source record has no center")
        elif not fpBbox.contains(center):
            raise Exception("The center is outside the Footprint of the source record")

        ccds = exposure.getInfo().getCoaddInputs().ccds
        measRecord.set(self.numberKey, len(ccds.subsetContaining(center, exposure.getWcs())))

    def fail(self, measRecord, error=None):
        measRecord.set(self.numberFlag, True)

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
        self.keyX = schema.addField(name + "_x", type="D", doc="peak centroid", units="pixels")
        self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixels")
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


class SingleFrameClassificationConfig(SingleFramePluginConfig):

    fluxRatio = lsst.pex.config.Field(dtype=float, default=.925, optional=True,
                                  doc="critical ratio of model to psf flux")
    modelErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                  doc="correction factor for modelFlux error")
    psfErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                  doc="correction factor for psfFlux error")

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

    @classmethod
    def getExecutionOrder(cls):
        return cls.CLASSIFY_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for extended sources, 0 for point sources.")
        self.keyFlag = schema.addField(name + "_flag", type="Flag", doc="Set to 1 for any fatal failure.")

    def measure(self, measRecord, exposure):
        modelFlux = measRecord.getModelFlux()
        psfFlux = measRecord.getPsfFlux()
        modelFluxFlag = (measRecord.getModelFluxFlag()
                         if measRecord.table.getModelFluxFlagKey().isValid()
                         else False)
        psfFluxFlag = (measRecord.getPsfFluxFlag()
                       if measRecord.table.getPsfFluxFlagKey().isValid()
                       else False)
        flux1 = self.config.fluxRatio*modelFlux
        if not self.config.modelErrFactor == 0:
            flux1 += self.config.modelErrFactor*measRecord.getModelFluxErr()
        flux2 = psfFlux
        if not self.config.psfErrFactor == 0:
            flux2 += self.config.psfErrFactor*measRecord.getPsfFluxErr()

        # A generic failure occurs when either FluxFlag is set to True
        # A generic failure also occurs if either calculated flux value is NAN:
        #     this can occur if the Flux field itself is NAN,
        #     or the ErrFactor != 0 and the FluxErr is NAN
        if numpy.isnan(flux1) or numpy.isnan(flux2) or modelFluxFlag or psfFluxFlag:
            self.fail(measRecord)
        else:
            if flux1 < flux2:
                measRecord.set(self.keyProbability, 0.0)
            else:
                measRecord.set(self.keyProbability, 1.0)

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception.  We should be setting a flag
        # instead.
        measRecord.set(self.keyFlag, True)


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
                                units="pixels^2")
        yyKey = schema.addField(name + "_yy", type="D", doc="transformed reference shape y^2 moment",
                                units="pixels^2")
        xyKey = schema.addField(name + "_xy", type="D", doc="transformed reference shape xy moment",
                                units="pixels^2")
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
