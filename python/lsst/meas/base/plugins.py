from lsst.meas.base.sfm import *
from lsst.meas.base.forced import *
from lsst.pex.exceptions.exceptionsLib import LsstCppException,LengthErrorException
from lsst.meas.base.baseLib import CentroidAlgorithmMapper,CentroidAlgorithmResult,FluxAlgorithmMapper,\
    FluxAlgorithmResult
import lsst.meas.base.base
import lsst.meas.base.baseLib
import lsst.afw.detection
import numpy
import lsst.meas.algorithms

class SingleFramePluginCpp(SingleFramePlugin):

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        any = lsst.meas.algorithms.AlgorithmRegistry.all.makeField(
            "field that accepts any algorithm"
        )
        find =  any.registry[name]
        control = find(find.ConfigClass())
        self.alg = control.makeAlgorithm(schema)

    def measureSingle(self, exposure, source):
        try:
            center = source.getCentroid()
            self.alg.apply(source, exposure, center)
        except LsstCppException as e:
            pass

    def measureMulti(self, exposure, sources):
        return


class ForcedPluginCpp(ForcedPlugin):

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        any = lsst.meas.algorithms.AlgorithmRegistry.all.makeField(
            "field that accepts any algorithm"
        )
        find =  any.registry[name]
        control = find(find.ConfigClass())
        schema = schemaMapper.getOutputSchema()
        beginlen = len(schema.asList())
        self.alg = control.makeAlgorithm(schema)
        newlist = schema.asList()
        for i in range(beginlen, len(newlist)):
            schemaMapper.addOutputField(newlist[i].getField())
        
    def measureSingle(self, exposure, source, refRecord, referenceWcs):
        try: 
            center = source.getCentroid()
            expWcs = exposure.getWcs()
            center = exposure.getWcs().skyToPixel(referenceWcs.pixelToSky(center))
            self.alg.apply(source, exposure, center)
        except Exception as e:
            pass #print "err in %s = %s" % (type(self), e.message)

    def measureMulti(self, exposure, sources):
        raise NotImplementedError()

class SingleFrameCentroidPlugin(SingleFramePlugin):
    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
 
class ForcedCentroidPlugin(ForcedPlugin):
    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.getSchema()
        beginlen = len(schema.asList)
        newlist = schema.asList()
        for i in range(beginlen, len(newlist)):
            schemaMapper.addOutputField(newlist[i].getField())
#--------------------------------------
# New sfm plugins
class SdssCentroidConfig(SingleFramePluginConfig):
    doMeasureSingle = True
    doMeasureMulti = False
    executionOrder = lsst.pex.config.Field(dtype=float, default=-1.0, doc="sets relative order of algorithms")

class SdssCentroid(SingleFramePluginCpp):

    ConfigClass = SdssCentroidConfig

    def measureSingle(self, exposure, source):
        try: 
            foot = source.getFootprint() 
            peak = foot.getPeaks()[0]
            center = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
            self.alg.apply(source, exposure, center)
        except LsstCppException as e:
            pass

SingleFramePlugin.registry.register("centroid.sdss", SdssCentroid)

class FluxGaussianConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")


class FluxGaussian(SingleFramePluginCpp):

    ConfigClass = FluxGaussianConfig

SingleFramePlugin.registry.register("flux.gaussian", FluxGaussian)


class FlagsPixelConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class FlagsPixel(SingleFramePluginCpp):

    ConfigClass = FlagsPixelConfig


SingleFramePlugin.registry.register("flags.pixel", FlagsPixel)

class CentroidGaussianConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class CentroidGaussian(SingleFramePluginCpp):

    ConfigClass = CentroidGaussianConfig


SingleFramePlugin.registry.register("centroid.gaussian", CentroidGaussian)

class CentroidNaiveConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class CentroidNaive(SingleFramePluginCpp):

    ConfigClass = CentroidNaiveConfig


SingleFramePlugin.registry.register("centroid.naive", CentroidNaive)


class FluxNaiveConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxNaive(SingleFramePluginCpp):

    ConfigClass = FluxNaiveConfig


SingleFramePlugin.registry.register("flux.naive", FluxNaive)

class FluxSincConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxSinc(SingleFramePluginCpp):

    ConfigClass = FluxSincConfig


SingleFramePlugin.registry.register("flux.sinc", FluxSinc)

class ClassificationExtendednessConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=5.0, doc="sets relative order of algorithms")

class ClassificationExtendedness(SingleFramePluginCpp):

    ConfigClass = ClassificationExtendednessConfig


SingleFramePlugin.registry.register("classification.extendedness", ClassificationExtendedness)

class SkyCoordConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=5.0, doc="sets relative order of algorithms")

class SkyCoord(SingleFramePluginCpp):

    ConfigClass = SkyCoordConfig

SingleFramePlugin.registry.register("skycoord", SkyCoord)

# ------------------------------------------------
# New forced plugins

class ForcedSdssCentroidConfig(ForcedPluginConfig):
    doMeasureSingle = True
    doMeasureMulti = False
    executionOrder = lsst.pex.config.Field(dtype=float, default=-1.0, doc="sets relative order of algorithms")

class ForcedSdssCentroid(ForcedPluginCpp):

    ConfigClass = ForcedSdssCentroidConfig

    def measureSingle(self, exposure, source, ref, refWcs):
        try: 
            foot = source.getFootprint() 
            peak = foot.getPeaks()[0]
            center = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
            self.alg.apply(source, exposure, center)
        except Exception as e:
            pass #print "err in sdss.centroid: %s, id = %d"%(e.message,source.getId())

ForcedPlugin.registry.register("centroid.sdss", ForcedSdssCentroid)

class ForcedFluxGaussianConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class ForcedFluxGaussian(ForcedPluginCpp):

    ConfigClass = ForcedFluxGaussianConfig


ForcedPlugin.registry.register("flux.gaussian", ForcedFluxGaussian)


class ForcedFlagsPixelConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class ForcedFlagsPixel(ForcedPluginCpp):

    ConfigClass = ForcedFlagsPixelConfig


ForcedPlugin.registry.register("flags.pixel", ForcedFlagsPixel)

class ForcedCentroidGaussianConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class ForcedCentroidGaussian(ForcedPluginCpp):

    ConfigClass = ForcedCentroidGaussianConfig


ForcedPlugin.registry.register("centroid.gaussian", ForcedCentroidGaussian)

class ForcedCentroidNaiveConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class ForcedCentroidNaive(ForcedPluginCpp):

    ConfigClass = ForcedCentroidNaiveConfig


ForcedPlugin.registry.register("centroid.naive", ForcedCentroidNaive)


class ForcedFluxNaiveConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class ForcedFluxNaive(ForcedPluginCpp):

    ConfigClass = ForcedFluxNaiveConfig


ForcedPlugin.registry.register("flux.naive", ForcedFluxNaive)

class ForcedFluxSincConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class ForcedFluxSinc(ForcedPluginCpp):

    ConfigClass = ForcedFluxSincConfig


ForcedPlugin.registry.register("flux.sinc", ForcedFluxSinc)

class ForcedClassificationExtendednessConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=5.0, doc="sets relative order of algorithms")

class ForcedClassificationExtendedness(ForcedPluginCpp):

    ConfigClass = ForcedClassificationExtendednessConfig


ForcedPlugin.registry.register("classification.extendedness", ForcedClassificationExtendedness)

class ForcedSkyCoordConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=5.0, doc="sets relative order of algorithms")

class ForcedSkyCoord(ForcedPluginCpp):

    ConfigClass = ForcedSkyCoordConfig

ForcedPlugin.registry.register("skycoord", ForcedSkyCoord)

# ----------------------------------------------
# New plugins
class SingleFramePeakCentroidConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class PeakCentroid(SingleFramePlugin):

    ConfigClass = SingleFramePeakCentroidConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        mapper = CentroidAlgorithmMapper(schema, name)
        self.config = config
        
    def measureSingle(self, exposure, source):
        peak = source.getFootprint().getPeaks()[0]
        result = lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())
        schema = source.getSchema()
        key = schema.find("centroid.peak").key
        source.set(key, result)

    def measureMulti(self, exposure, sources):
        raise NotImplementedError()


SingleFramePlugin.registry.register("centroid.peak", PeakCentroid)

class FluxPsfConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class FluxPsf(SingleFramePlugin):

    ConfigClass = FluxPsfConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.config = config
        self.algorithm = lsst.meas.base.baseLib.PsfFluxAlgorithm()
        self.resultMapper = self.algorithm.makeResultMapper(schema)
        
    def measureSingle(self, exposure, source):
        try:
            center = source.getCentroid()
            result = self.algorithm.apply(exposure, source.getFootprint(), center)
            self.resultMapper.apply(source, result)

        except LsstCppException as e:
            pass

    def measureMulti(self, exposure, sources):
        raise NotImplementedError()


SingleFramePlugin.registry.register("flux.psf", FluxPsf)


class SdssShapeConfig(SingleFramePluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class SdssShape(SingleFramePlugin):

    ConfigClass = FluxPsfConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.config = config
        self.control = lsst.meas.base.baseLib.SdssShapeControl()
        self.algorithm = lsst.meas.base.baseLib.SdssShapeAlgorithm()
        self.resultMapper = self.algorithm.makeResultMapper(schema)
        
    def measureSingle(self, exposure, source):
        try: 
            center = source.getCentroid()
            result = self.algorithm.apply(self.control, exposure.getMaskedImage(),
                source.getFootprint(), center)
            self.resultMapper.apply(source, result)

        except LsstCppException as e:
            pass
        except Exception as e:
            print "Raise bad exception"
            raise e
            
    def measureMulti(self, exposure, sources):
        raise NotImplementedError()


SingleFramePlugin.registry.register("shape.sdss", SdssShape)

# ----------------------------------------------

class ForcedCentroidPlugin(ForcedPlugin):
    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        schema = schemaMapper.getOutputSchema()
        beginlen = len(schema.asList())
        CentroidAlgorithmMapper(schema, name)   # add the names of the output fields to the schema
        newlist = schema.asList()
        for i in range(beginlen, len(newlist)):
            schemaMapper.addOutputField(newlist[i].getField())

class ForcedPeakCentroidConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=0.0, doc="sets relative order of algorithms")

class PeakCentroid(ForcedCentroidPlugin):

    ConfigClass = ForcedPeakCentroidConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedCentroidPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        self.name = name
        self.config = config
        
    def measureSingle(self, exposure, source, refRecord, referenceWcs):
        targetWcs = exposure.getWcs()       
        peak = source.getFootprint().getPeaks()[0]
        result = targetWcs.skyToPixel(
            referenceWcs.pixelToSky(lsst.afw.geom.Point2D(peak.getFx(), peak.getFy())))
        key = source.getSchema().find("centroid.peak").key
        source.set(key, result)
        
    def measureMulti(self, exposure, sources):
        raise NotImplementedError()


ForcedPlugin.registry.register("centroid.peak", PeakCentroid)

class ForcedPsfFluxConfig(ForcedPluginConfig):
    executionOrder = lsst.pex.config.Field(dtype=float, default=2.0, doc="sets relative order of algorithms")

class ForcedPsfFlux(ForcedPlugin):

    ConfigClass = ForcedPsfFluxConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        self.algorithm = lsst.meas.base.baseLib.PsfFluxAlgorithm()
        schema = schemaMapper.getOutputSchema()
        self.resultMapper = self.algorithm.makeResultMapper(schema)
        beginlen = len(schema.asList())
        FluxAlgorithmMapper(schema, name)   # add the names of the output fields to the schema
        newlist = schema.asList()
        for i in range(beginlen, len(newlist)):
            schemaMapper.addOutputField(newlist[i].getField())
        
    def measureSingle(self, exposure, source, refRecord, referenceWcs):
        try:
            center = source.getCentroid()
            expregion = lsst.afw.geom.Box2I(exposure.getXY0(),
                lsst.afw.geom.Extent2I(exposure.getWidth(), exposure.getHeight()))
            footprint = source.getFootprint()
            if footprint.isHeavy(): footprint = lsst.afw.detection.Footprint(footprint)
            result = self.algorithm.apply(exposure, 
                footprint.transform(referenceWcs, exposure.getWcs(), expregion, True), center)
            self.resultMapper.apply(source, result)
        except LsstCppException as e:
            pass

    def measureMulti(self, exposure, sources):
        raise NotImplementedError()

ForcedPlugin.registry.register("flux.psf", ForcedPsfFlux)

class ForcedSdssShapeConfig(ForcedPluginConfig):
    doMeasureSingle = True
    doMeasureMulti = False
    executionOrder = lsst.pex.config.Field(dtype=float, default=1.0, doc="sets relative order of algorithms")



class ForcedSdssShape(ForcedPlugin):

    ConfigClass = ForcedSdssShapeConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        self.algorithm = lsst.meas.base.baseLib.SdssShapeAlgorithm()
        schema = schemaMapper.getOutputSchema()
        self.resultMapper = self.algorithm.makeResultMapper(schema)
        beginlen = len(schema.asList())
        SdssShapeAlgorithmMapper(schema, name)   # add the names of the output fields to the schema
        newlist = schema.asList()
        for i in range(beginlen, len(newlist)):
            schemaMapper.addOutputField(newlist[i].getField())
        
    def measureSingle(self, exposure, source, refRecord, referenceWcs):
        try:
            center = source.getCentroid()
            expregion = lsst.afw.geom.Box2I(exposure.getXY0(),
                lsst.afw.geom.Extent2I(exposure.getWidth(), exposure.getHeight()))
            footprint = source.getFootprint()
            if footprint.isHeavy(): footprint = lsst.afw.detection.Footprint(footprint)
            footprint = footprint.transform(referenceWcs, exposure.getWcs(), expregion, True)
            result = self.algorithm.apply(self.control, exposure.getMaskedImage(),
                footprint, center)
            self.resultMapper.apply(source, result)
        except LsstCppException as e:
            pass

    def measureMulti(self, exposure, sources):
        raise NotImplementedError()


ForcedPlugin.registry.register("shape.sdss", ForcedSdssShape)
