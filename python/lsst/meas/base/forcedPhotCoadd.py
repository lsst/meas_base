# This file is part of meas_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import lsst.pex.config
import lsst.pipe.base
import lsst.coadd.utils
import lsst.afw.table

from .forcedPhotImage import ForcedPhotImageConfig, ForcedPhotImageTask

__all__ = ("ForcedPhotCoaddConfig", "ForcedPhotCoaddTask")


class ForcedPhotCoaddConfig(ForcedPhotImageConfig):
    footprintDatasetName = lsst.pex.config.Field(
        doc="Dataset (without coadd prefix) that should be used to obtain (Heavy)Footprints for sources. "
            "Must have IDs that match those of the reference catalog."
            "If None, Footprints will be generated by transforming the reference Footprints.",
        dtype=str,
        default="meas",
        optional=True
    )

    hasFakes = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Should be set to True if fake sources have been inserted into the input data."
    )

    def setDefaults(self):
        ForcedPhotImageTask.ConfigClass.setDefaults(self)
        # Copy 'id' and 'parent' columns without renaming them; since these are
        # the only IDs we have in coadd processing, there's no need to qualify
        # them with 'object'.
        self.measurement.copyColumns["id"] = "id"
        self.measurement.copyColumns["parent"] = "parent"
        self.references.removePatchOverlaps = False  # see validate() for why
        self.measurement.plugins.names |= ['base_InputCount', 'base_Variance']
        self.measurement.plugins['base_PixelFlags'].masksFpAnywhere = ['CLIPPED', 'SENSOR_EDGE',
                                                                       'REJECTED', 'INEXACT_PSF']
        self.measurement.plugins['base_PixelFlags'].masksFpCenter = ['CLIPPED', 'SENSOR_EDGE',
                                                                     'REJECTED', 'INEXACT_PSF']

    def validate(self):
        ForcedPhotImageTask.ConfigClass.validate(self)
        if (self.measurement.doReplaceWithNoise and self.footprintDatasetName is not None and
                self.references.removePatchOverlaps):
            raise ValueError("Cannot use removePatchOverlaps=True with deblended footprints, as parent "
                             "sources may be rejected while their children are not.")


class ForcedPhotCoaddRunner(lsst.pipe.base.ButlerInitializedTaskRunner):
    """Get the psfCache setting into ForcedPhotCoaddTask"""
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return lsst.pipe.base.ButlerInitializedTaskRunner.getTargetList(parsedCmd,
                                                                        psfCache=parsedCmd.psfCache)


class ForcedPhotCoaddTask(ForcedPhotImageTask):
    """A command-line driver for performing forced measurement on coadd images.

    Notes
    -----
    In addition to the run method, `ForcedPhotCcdTask` overrides several
    methods of `ForcedPhotImageTask` to specialize it for coadd processing,
    including `~ForcedPhotImageTask.makeIdFactory` and
    `~ForcedPhotImageTask.fetchReferences`. None of these should be called
    directly by the user, though it may be useful to override them further in
    subclasses.
    """

    ConfigClass = ForcedPhotCoaddConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCoadd"
    dataPrefix = "deepCoadd_"

    def getExposure(self, dataRef):

        if self.config.hasFakes:
            name = "fakes_" + self.config.coaddName + "Coadd_calexp"
        else:
            name = self.config.coaddName + "Coadd_calexp"

        return dataRef.get(name) if dataRef.datasetExists(name) else None

    def makeIdFactory(self, dataRef):
        """Create an object that generates globally unique source IDs.

        Source IDs are created based on a per-CCD ID and the ID of the CCD
        itself.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference. The "CoaddId_bits" and "CoaddId" datasets
            are accessed. The data ID must have tract and patch keys.
        """
        # With the default configuration, this IdFactory doesn't do anything,
        # because the IDs it generates are immediately overwritten by the ID
        # from the reference catalog (since that's in
        # config.measurement.copyColumns).  But we create one here anyway, to
        # allow us to revert back to the old behavior of generating new forced
        # source IDs, just by renaming the ID in config.copyColumns to
        # "object_id".
        expBits = dataRef.get(self.config.coaddName + "CoaddId_bits")
        expId = int(dataRef.get(self.config.coaddName + "CoaddId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    def getExposureId(self, dataRef):
        return int(dataRef.get(self.config.coaddName + "CoaddId"))

    def fetchReferences(self, dataRef, exposure):
        """Return an iterable of reference sources which overlap the exposure.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference corresponding to the image to be measured;
            should have tract, patch, and filter keys.

        exposure : `lsst.afw.image.Exposure`
            Unused.

        Notes
        -----
        All work is delegated to the references subtask; see
        `CoaddSrcReferencesTask` for information about the default behavior.
        """
        skyMap = dataRef.get(self.dataPrefix + "skyMap", immediate=True)
        tractInfo = skyMap[dataRef.dataId["tract"]]
        patch = tuple(int(v) for v in dataRef.dataId["patch"].split(","))
        patchInfo = tractInfo.getPatchInfo(patch)
        references = lsst.afw.table.SourceCatalog(self.references.schema)
        references.extend(self.references.fetchInPatches(dataRef, patchList=[patchInfo]))
        return references

    def attachFootprints(self, sources, refCat, exposure, refWcs, dataRef):
        r"""Attach Footprints to source records.

        For coadd forced photometry, we use the deblended "heavy"
        `~lsst.afw.detection.Footprint`\ s from the single-band measurements
        of the same band - because we've guaranteed that the peaks (and hence
        child sources) will be consistent across all bands before we get to
        measurement, this should yield reasonable deblending for most sources.
        It's most likely limitation is that it will not provide good flux
        upper limits for sources that were not detected in this band but were
        blended with sources that were.
        """
        if self.config.footprintDatasetName is None:
            return ForcedPhotImageTask.attachFootprints(self, sources, refCat, exposure, refWcs, dataRef)
        self.log.info("Loading deblended footprints for sources from %s, %s" %
                      (self.config.footprintDatasetName, dataRef.dataId))
        fpCat = dataRef.get("%sCoadd_%s" % (self.config.coaddName, self.config.footprintDatasetName),
                            immediate=True)
        for refRecord, srcRecord in zip(refCat, sources):
            fpRecord = fpCat.find(refRecord.getId())
            if fpRecord is None:
                raise LookupError("Cannot find Footprint for source %s; please check that %sCoadd_%s "
                                  "IDs are compatible with reference source IDs" %
                                  (srcRecord.getId(), self.config.coaddName,
                                   self.config.footprintDatasetName))
            srcRecord.setFootprint(fpRecord.getFootprint())

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=lsst.coadd.utils.CoaddDataIdContainer)
        parser.add_argument("--psfCache", type=int, default=100, help="Size of CoaddPsf cache")
        return parser
