"""Class for force measuring a Coadd exposure using a reference catalog made from measuring Coadds.
   See comments in the parent class in forced.py for more information.
"""

import lsst.pex.config
from lsst.pipe.base import Task, CmdLineTask, Struct, timeMethod, ArgumentParser, ButlerInitializedTaskRunner
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer
import lsst.daf.base
from lsst.pex.config import DictField,ConfigurableField
from .forced import *
from .base import *
from .references import CoaddSrcReferencesTask
__all__ = ("ForcedCoaddMeasurementConfig", "ForcedCoaddMeasurementTask")

class ForcedCoaddMeasurementConfig(ForcedMeasurementConfig):
    pass


class ForcedCoaddMeasurementTask(ForcedMeasurementTask):
    """Forced measurement driver task

    This task is intended as a command-line script base class, in the model of ProcessImageTask
    (i.e. it should be subclasses for running on CCDs and Coadds).
    """

    ConfigClass = ForcedCoaddMeasurementConfig
    RunnerClass = ButlerInitializedTaskRunner
    _DefaultName = "forcedCoaddMeasurementTask"
    dataPrefix = "deepCoadd_"

    # Find references which overlap a skyMap made up of one or more patches.
    # The references in each patch correspond to measurements from individual coadss.       
    def fetchReferences(self, dataRef, exposure):
        skyMap = dataRef.get(self.dataPrefix + "skyMap", immediate=True)
        tractInfo = skyMap[dataRef.dataId["tract"]]
        patch = tuple(int(v) for v in dataRef.dataId["patch"].split(","))
        patchInfo = tractInfo.getPatchInfo(patch)
        return self.references.fetchInPatches(dataRef, patchList=[patchInfo])

    def makeIdFactory(self, dataRef):
        expBits = dataRef.get(self.config.coaddName + "CoaddId_bits")
        expId = long(dataRef.get(self.config.coaddName + "CoaddId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    def writeOutput(self, dataRef, sources):
        """Write forced source table

        @param dataRef  Data reference from butler
        @param sources  SourceCatalog to save
        """
        dataRef.put(sources, self.dataPrefix + "forced_src")

    def getExposure(self, dataRef):
        """Read input exposure on which to perform the measurements

        @param dataRef       Data reference from butler
        """
        return dataRef.get(self.dataPrefix + "calexp", immediate=True)

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def getPreviousTaskClass(self):
        return MeasureCoaddTask

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%sCoadd_forced_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%sCoadd_forced_metadata" % (self.config.coaddName,)
