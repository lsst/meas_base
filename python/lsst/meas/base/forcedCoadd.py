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
   See comments in the parent class in forcedImage.py for more information.
"""

import lsst.pex.config

from lsst.pipe.base import Task, CmdLineTask, Struct, timeMethod, ArgumentParser, ButlerInitializedTaskRunner
from lsst.coadd.utils import CoaddDataIdContainer
import lsst.daf.base
from lsst.pex.config import DictField,ConfigurableField

from .forcedImage import *
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
        """An Id Factory is used to convert ids from the previous pipelines into an id,
        using specific bit field which the dataRef butler knows about (camera specific)
        """
        expBits = dataRef.get(self.config.coaddName + "CoaddId_bits")
        expId = long(dataRef.get(self.config.coaddName + "CoaddId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)


    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=CoaddDataIdContainer)
        return parser


