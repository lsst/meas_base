"""Base classes for forced measurement plugin algorithms and the driver task for these.

In forced measurement, a reference catalog is used to define restricted measurements (usually just fluxes)
on an image.  As the reference catalog may be deeper than the detection limit of the measurement image, we
do not assume that we can use detection and deblend information from the measurement image.  Instead, we
assume this information is present in the reference catalog and has been "transformed" in some sense to
the measurement frame.  At the very least, this means that Footprints from the reference catalog should
be transformed and installed as Footprints in the output measurement catalog.  If we have a procedure that
can transform HeavyFootprints, we can then proceed with measurement as usual, but using the reference
catalog's id and parent fields to define deblend families.  If this transformation does not preserve
HeavyFootprints (this is currently the case), then we will only be able to replace objects with noise
one deblend family at a time, and hence measurements run in single-object mode may be contaminated by
neighbors when run on objects with parent != 0.

Measurements are generally recorded in the coordinate system of the image being measured (and all
slot-eligible fields must be), but non-slot fields may be recorded in other coordinate systems if necessary
to avoid information loss (this should, of course, be indicated in the field documentation).  Note that
the reference catalog may be in a different coordinate system; it is the responsibility of algorithms
to transform the data they need themselves, using the reference WCS provided.  However, for algorithms
that only require a position, they may simply use output SourceCatalog's centroid slot, which will generally
be set to the transformed position of the reference object before any algorithms are run, and hence avoid
using the reference catalog at all.
"""

#
# NOTE: to fully implement forced measurement, we'll need not just this measurement task, but also a
# command-line driver task that uses it as a subtask.  A good prototype for this already exists on
# the HSC branch of pipe_tasks, in the forcedPhot*.py and references.py files.  These should be
# transferred to the LSST side of pipe_tasks, and modified to use the new measurement API described
# here.  Note that it will also be responsible for transforming Footprints (in the generateSources())
# method, and attaching them to sources, as this is no longer something measurement tasks do.  There
# are also some hackish workarounds in that code that could be removed with properly vetted
# modifications to pipe_base (i.e. we need to get a Butler in Task ctors).  We should do that now
# as well.
#
# It's also worth considering merging that command-line-driver task with this measurement task, and
# flattening the hierarchy.  In that case, though, we'd probably want to put the base command-line
# task in meas_base and let the CCD- and coadd-specialized subclasses continue to be in pipe_tasks
# to avoid dependency issues.
#

import lsst.pex.config
from lsst.pipe.base import Task, CmdLineTask, Struct, timeMethod, ArgumentParser, ButlerInitializedTaskRunner
import lsst.daf.base
from lsst.pex.config import DictField,ConfigurableField
from lsst.pipe.tasks.references import CoaddSrcReferencesTask
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer
from .forced import *
from .base import *
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


