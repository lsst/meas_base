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

"""Class for force measuring a Coadd exposure using a reference catalog made from measuring Coadds.
"""

import lsst.pex.config
import lsst.pipe.base
import lsst.coadd.utils
import lsst.afw.table

from .forcedImage import *
from .base import *

__all__ = ("ProcessForcedCoaddConfig", "ProcessForcedCoaddTask")

class ProcessForcedCoaddConfig(ProcessImageForcedConfig):
    pass

## \addtogroup LSST_task_documentation
## \{
## \page processForcedCoaddTask
## \ref ProcessForcedCoaddTask_ "ProcessForcedCoaddTask"
## \copybrief ProcessForcedCoaddTask
## \}

class ProcessForcedCoaddTask(ProcessImageForcedTask):
    """!
    \anchor ProcessForcedCoaddTask_

    \brief The ProcessForcedCoaddTask is used to measure the properties of sources on a Coadd, using
           sources from a run on another coadd (usually a different band) as references.

    \section meas_base_processForcedCoaddTask_Contents Contents

     - \ref meas_base_processForcedCoaddTask_Purpose
     - \ref meas_base_processForcedCoaddTask_IO
     - \ref meas_base_processForcedCoaddTask_Config

    \section meas_base_processForcedCoaddTask_Purpose	Description

    This task is a subclass of ProcessForcedImageTask which is specifically for doing forced
    measurement on a coadd, using the as a reference catalog the detections which
    were made on overlapping coadds.

    \section meas_base_processForcedCoaddTask_IO		Inputs/Outputs to the run method

    See \link lsst.meas.base.forcedImage.ProcessImageForcedTask#run ProcessImageForcedTask.run\endlink
    for more information.

    The run method (inherited from ProcessForcedImageTask) takes a lsst.daf.persistence.ButlerDataRef
    argument that corresponds to a coadd image.  This is used to provide all the inputs and outputs
    for the task:
     - A "*Coadd_src" (e.g. "deepCoadd_src") dataset is used as the reference catalog.  This not loaded
       directly from the passed dataRef, however; only the patch and tract are used, while the filter
       is set by the configuration for the references subtask (see CoaddSrcReferencesTask).
     - A "*Coadd_calexp" (e.g. "deepCoadd_calexp") dataset is used as the measurement image.  Note that
       this means that ProcessCoaddTask must be run on an image before ProcessForcedCoaddTask, in order
       to generate the "*Coadd_calexp" dataset.
     - A "*Coadd_forced_src" (e.g. "deepCoadd_forced_src") dataset will be written with the output
       measurement catalog.

    In addition to the run method, ProcessForcedCcdTask overrides several methods of ProcessForcedImageTask
    to specialize it for coadd processing, including makeIdFactory() and fetchReferences().  None of these
    should be called directly by the user, though it may be useful to override them further in subclasses.

    \section meas_base_processForcedCoaddTask_Config       Configuration parameters

    See \ref ProcessForcedCoaddConfig
    """

    ConfigClass = ProcessForcedCoaddConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedCoaddTask"
    dataPrefix = "deepCoadd_"

    def fetchReferences(self, dataRef, exposure):
        """Find references which overlap a skyMap made up of one or more patches.
        The references in each patch correspond to measurements from individual coadss.
        """
        skyMap = dataRef.get(self.dataPrefix + "skyMap", immediate=True)
        tractInfo = skyMap[dataRef.dataId["tract"]]
        patch = tuple(int(v) for v in dataRef.dataId["patch"].split(","))
        patchInfo = tractInfo.getPatchInfo(patch)
        return self.references.fetchInPatches(dataRef, patchList=[patchInfo])

    def makeIdFactory(self, dataRef):
        """An IdFactory is used to convert ids from the previous pipelines into an id,
        using specific bit field which the dataRef butler knows about (camera specific)
        """
        expBits = dataRef.get(self.config.coaddName + "CoaddId_bits")
        expId = long(dataRef.get(self.config.coaddName + "CoaddId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=lsst.coadd.utils.CoaddDataIdContainer)
        return parser
