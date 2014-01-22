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


    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=CoaddDataIdContainer)
        return parser


