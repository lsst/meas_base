#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2014 LSST Corporation.
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

import lsst.pex.config
import lsst.pipe.base
import lsst.coadd.utils
import lsst.afw.table

from .forcedPhotImage import *
from .base import *
from .references import CoaddSrcReferencesTask

__all__ = ("ForcedPhotCoaddConfig", "ForcedPhotCoaddTask")

class ForcedPhotCoaddConfig(ProcessImageForcedConfig):
    pass

## @addtogroup LSST_task_documentation
## @{
## @page processForcedCoaddTask
## ForcedPhotCoaddTask
## @copybrief ForcedPhotCoaddTask
## @}

class ForcedPhotCoaddTask(ProcessImageForcedTask):
    """!
    A command-line driver for performing forced measurement on coadd images

    This task is a subclass of ForcedPhotImageTask which is specifically for doing forced
    measurement on a coadd, using as a reference catalog detections which were made on overlapping
    coadds (i.e. in other bands).

    The run method (inherited from ForcedPhotImageTask) takes a lsst.daf.persistence.ButlerDataRef
    argument that corresponds to a coadd image.  This is used to provide all the inputs and outputs
    for the task:
     - A "*Coadd_src" (e.g. "deepCoadd_src") dataset is used as the reference catalog.  This not loaded
       directly from the passed dataRef, however; only the patch and tract are used, while the filter
       is set by the configuration for the references subtask (see CoaddSrcReferencesTask).
     - A "*Coadd_calexp" (e.g. "deepCoadd_calexp") dataset is used as the measurement image.  Note that
       this means that ProcessCoaddTask must be run on an image before ForcedPhotCoaddTask, in order
       to generate the "*Coadd_calexp" dataset.
     - A "*Coadd_forced_src" (e.g. "deepCoadd_forced_src") dataset will be written with the output
       measurement catalog.

    In addition to the run method, ForcedPhotCcdTask overrides several methods of ForcedPhotImageTask
    to specialize it for coadd processing, including makeIdFactory() and fetchReferences().  None of these
    should be called directly by the user, though it may be useful to override them further in subclasses.
    """

    ConfigClass = ForcedPhotCoaddConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCoaddTask"
    dataPrefix = "deepCoadd_"

    def makeIdFactory(self, dataRef):
        """Create an object that generates globally unique source IDs from per-CCD IDs and the CCD ID.

        @param dataRef       Data reference from butler.  The "CoaddId_bits" and "CoaddId"
                             datasets are accessed.  The data ID must have tract and patch keys.
        """
        expBits = dataRef.get(self.config.coaddName + "CoaddId_bits")
        expId = long(dataRef.get(self.config.coaddName + "CoaddId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    def fetchReferences(self, dataRef, exposure):
        """Return an iterable of reference sources which overlap the exposure

        @param dataRef       Data reference from butler corresponding to the image to be measured;
                             should have tract, patch, and filter keys.
        @param exposure      lsst.afw.image.Exposure to be measured (not used by this implementation)

        All work is delegated to the references subtask; see CoaddSrcReferencesTask for information
        about the default behavior.
        """
        skyMap = dataRef.get(self.dataPrefix + "skyMap", immediate=True)
        tractInfo = skyMap[dataRef.dataId["tract"]]
        patch = tuple(int(v) for v in dataRef.dataId["patch"].split(","))
        patchInfo = tractInfo.getPatchInfo(patch)
        return self.references.fetchInPatches(dataRef, patchList=[patchInfo])

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=lsst.coadd.utils.CoaddDataIdContainer)
        return parser


