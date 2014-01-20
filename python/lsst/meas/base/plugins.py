from lsst.meas.base.sfm import *
from lsst.meas.base.forced import *
from lsst.afw.table import *
import lsst.afw.table.tableLib as tableLib
from lsst.pex.exceptions.exceptionsLib import LsstCppException,LengthErrorException
import lsst.meas.base.base
import lsst.afw.detection
import numpy
import lsst.meas.algorithms

# ----------------------------------------------
# New plugins written purely in Python, based on meas_base plugin model

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
        schema.addField("centroid.peak.cov", "CovPointF", doc="covariance of measured centroid", units="pixels^2")        
    def measureSingle(self, exposure, source):
        peak = source.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        schema = source.getSchema()
        key = schema.find("centroid.peak").key
        source.set(key, result)

    def measureMulti(self, exposure, sources):
        raise NotImplementedError()

SingleFramePlugin.registry.register("centroid.peak", PeakCentroid)
# ----------------------------------------------

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
        
    def measureSingle(self, exposure, source, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()       
        peak = refRecord.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not referenceWcs == targetWcs: 
            result = targetWcs.skyToPixel(referenceWcs.pixelToSky(result))
        source.set(self.peakkey, result)
        
    def measureMulti(self, exposure, sources):
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
        field = tableLib.Field_CovPointF("transformed.peak.cov", "measured peak centroid covariance", "pixels")
        self.covkey = schemaMapper.addOutputField(field)
        
    def measureSingle(self, exposure, source, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()      
        footprint = refRecord.getFootprint() 
        peak = footprint.getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        if not targetWcs == referenceWcs:
            result = targetWcs.skyToPixel(referenceWcs.pixelToSky(result))
        source.set(self.peakkey, result)
        
    def measureMulti(self, exposure, sources):
        raise NotImplementedError()


ForcedPlugin.registry.register("centroid.transformed", ForcedTransformedCentroid)

