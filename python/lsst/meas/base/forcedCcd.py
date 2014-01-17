"""Class to measure a Ccd calexp exposure, using a reference catalog from a Coadd
   See the comment in the parent class in forced.py for more information.
"""

import argparse
import lsst.pex.config
from lsst.pipe.base import Task, CmdLineTask, Struct, timeMethod, ArgumentParser, ButlerInitializedTaskRunner
import lsst.daf.base
from lsst.pex.config import DictField,ConfigurableField
from .forced import *
from .base import *

import lsst.afw.table
import lsst.afw.image
import lsst.pipe.base
from lsst.pex.config import Config, ConfigurableField, DictField, Field

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("ForcedCcdMeasurementConfig", "ForcedCcdMeasurementTask")

class ForcedCcdDataIdContainer(lsst.pipe.base.DataIdContainer):
    """A version of lsst.pipe.base.DataIdContainer specialized for forced photometry on CCDs.
    
    Required because we need to add "tract" to the raw data ID keys, and that's tricky.
    This IdContainer assumes that a calexp is being measured using the detection information
    from the set of coadds which intersect with the calexp.  a set of reference catalog
    from a coadd which overlaps it.  It needs the calexp id (visit, raft, sensor), but it
    also uses the tract to decide what set of coadds to use.  The references from the tract
    whose patches intersect with the calexp are used.
    """
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        for dataId in self.idList:
            if "tract" not in dataId:
                raise argparse.ArgumentError(None, "--id must include tract")
            tract = dataId.pop("tract")
            # making a DataRef for src fills out any missing keys and allows us to iterate
            for srcDataRef in namespace.butler.subset("src", dataId=dataId):
                forcedDataId = srcDataRef.dataId.copy()
                forcedDataId['tract'] = tract
                dataRef = namespace.butler.dataRef(
                    datasetType = "forced_src",
                    dataId = forcedDataId,
                    )
                self.refList.append(dataRef)

class ForcedCcdMeasurementConfig(ForcedMeasurementConfig):
    doApplyUberCal = Field(
        dtype = bool,
        doc = "Apply meas_mosaic ubercal results to input calexps?",
        default = True
    )


class ForcedCcdMeasurementTask(ForcedMeasurementTask):
    """Forced measurement driver task

    This task is intended as a command-line script base class, in the model of ProcessImageTask
    (i.e. it should be subclasses for running on Coadds and Ccds).
    """

    ConfigClass = ForcedCcdMeasurementConfig
    RunnerClass = ButlerInitializedTaskRunner
    _DefaultName = "forcedCcdMeasurementTask"
    dataPrefix = ""

    # For the output table, make a source id for each output from the 
    # exposure and ccd ids  
    def makeIdFactory(self, dataRef):
        expBits = dataRef.get("ccdExposureId_bits")
        expId = long(dataRef.get("ccdExposureId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)        

    # get the references which overlap the exposure
    def fetchReferences(self, dataRef, exposure):
        bbox = exposure.getBBox(lsst.afw.image.PARENT)
        return self.references.fetchInBox(dataRef, bbox, exposure.getWcs())

    def getExposure(self, dataRef):
        """Read input exposure to measure

        @param dataRef       Data reference from butler
        """
        exposure = dataRef.get("calexp", immediate=True)
        if not self.config.doApplyUberCal:
            return exposure
        if applyMosaicResults is None:
            raise RuntimeError(
                "Cannot use improved calibrations for %s because meas_mosaic could not be imported."
                % dataRef.dataId
                )
        else:
            applyMosaicResults(dataRef, calexp=exposure)
        return exposure

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=ForcedCcdDataIdContainer)
        return parser

