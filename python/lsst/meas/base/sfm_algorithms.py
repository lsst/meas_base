from lsst.meas.base.sfm import *
import lsst.afw.detection
import numpy

class TestCentroidConfig(SingleFrameAlgorithmConfig):
    fractional = lsst.pex.config.Field(dtype=bool, default=True,
                    doc="whether to center on fractional pixels")


#  Algoritm which does nothing
class TestNull(SingleFrameAlgorithm):

    ConfigClass = SingleFrameAlgorithmConfig
    doMeasureSingle = True
    doMeasureMulti = False
    def __init__(self, config, name, schema=None, flags=None, others=None, metadata=None):
        pass

    def measureSingle(self, exposure, source):
        return

    def measureMulti(self, exposure, sources):
        return

class TestFluxConfig(SingleFrameAlgorithmConfig):
    pass

SingleFrameAlgorithm.registry.register("test.null", TestNull)
