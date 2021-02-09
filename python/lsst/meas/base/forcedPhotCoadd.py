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
import lsst.coadd.utils
import lsst.afw.table

import lsst.pipe.base as pipeBase

from .references import MultiBandReferencesTask
from .forcedMeasurement import ForcedMeasurementTask
from .applyApCorr import ApplyApCorrTask
from .catalogCalculation import CatalogCalculationTask

__all__ = ("ForcedPhotCoaddConfig", "ForcedPhotCoaddTask")


class ForcedPhotCoaddRunner(pipeBase.ButlerInitializedTaskRunner):
    """Get the psfCache setting into ForcedPhotCoaddTask"""
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return pipeBase.ButlerInitializedTaskRunner.getTargetList(parsedCmd,
                                                                  psfCache=parsedCmd.psfCache)


class ForcedPhotCoaddConnections(pipeBase.PipelineTaskConnections,
                                 dimensions=("band", "skymap", "tract", "patch"),
                                 defaultTemplates={"inputCoaddName": "deep",
                                                   "outputCoaddName": "deep"}):
    inputSchema = pipeBase.connectionTypes.InitInput(
        doc="Schema for the input measurement catalogs.",
        name="{inputCoaddName}Coadd_ref_schema",
        storageClass="SourceCatalog",
    )
    outputSchema = pipeBase.connectionTypes.InitOutput(
        doc="Schema for the output forced measurement catalogs.",
        name="{outputCoaddName}Coadd_forced_src_schema",
        storageClass="SourceCatalog",
    )
    exposure = pipeBase.connectionTypes.Input(
        doc="Input exposure to perform photometry on.",
        name="{inputCoaddName}Coadd",
        storageClass="ExposureF",
        dimensions=["band", "skymap", "tract", "patch"],
    )
    refCat = pipeBase.connectionTypes.Input(
        doc="Catalog of shapes and positions at which to force photometry.",
        name="{inputCoaddName}Coadd_ref",
        storageClass="SourceCatalog",
        dimensions=["skymap", "tract", "patch"],
    )
    refCatInBand = pipeBase.connectionTypes.Input(
        doc="Catalog of shapes and positions in the band having forced photometry done",
        name="{inputCoaddName}Coadd_meas",
        storageClass="SourceCatalog",
        dimensions=("band", "skymap", "tract", "patch")
    )
    refWcs = pipeBase.connectionTypes.Input(
        doc="Reference world coordinate system.",
        name="{inputCoaddName}Coadd.wcs",
        storageClass="Wcs",
        dimensions=["band", "skymap", "tract", "patch"],
    )  # used in place of a skymap wcs because of DM-28880
    measCat = pipeBase.connectionTypes.Output(
        doc="Output forced photometry catalog.",
        name="{outputCoaddName}Coadd_forced_src",
        storageClass="SourceCatalog",
        dimensions=["band", "skymap", "tract", "patch"],
    )


class ForcedPhotCoaddConfig(pipeBase.PipelineTaskConfig,
                            pipelineConnections=ForcedPhotCoaddConnections):
    references = lsst.pex.config.ConfigurableField(
        target=MultiBandReferencesTask,
        doc="subtask to retrieve reference source catalog"
    )
    measurement = lsst.pex.config.ConfigurableField(
        target=ForcedMeasurementTask,
        doc="subtask to do forced measurement"
    )
    coaddName = lsst.pex.config.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )
    doApCorr = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Run subtask to apply aperture corrections"
    )
    applyApCorr = lsst.pex.config.ConfigurableField(
        target=ApplyApCorrTask,
        doc="Subtask to apply aperture corrections"
    )
    catalogCalculation = lsst.pex.config.ConfigurableField(
        target=CatalogCalculationTask,
        doc="Subtask to run catalogCalculation plugins on catalog"
    )
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
        # Docstring inherited.
        # Make catalogCalculation a no-op by default as no modelFlux is setup by default in
        # ForcedMeasurementTask
        super().setDefaults()

        self.catalogCalculation.plugins.names = []
        self.measurement.copyColumns["id"] = "id"
        self.measurement.copyColumns["parent"] = "parent"
        self.references.removePatchOverlaps = False  # see validate() for why
        self.measurement.plugins.names |= ['base_InputCount', 'base_Variance']
        self.measurement.plugins['base_PixelFlags'].masksFpAnywhere = ['CLIPPED', 'SENSOR_EDGE',
                                                                       'REJECTED', 'INEXACT_PSF']
        self.measurement.plugins['base_PixelFlags'].masksFpCenter = ['CLIPPED', 'SENSOR_EDGE',
                                                                     'REJECTED', 'INEXACT_PSF']

    def validate(self):
        super().validate()
        if (self.measurement.doReplaceWithNoise and self.footprintDatasetName is not None
                and self.references.removePatchOverlaps):
            raise ValueError("Cannot use removePatchOverlaps=True with deblended footprints, as parent "
                             "sources may be rejected while their children are not.")


class ForcedPhotCoaddTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    """A command-line driver for performing forced measurement on coadd images.

    Parameters
    ----------
    butler : `lsst.daf.persistence.butler.Butler`, optional
        A Butler which will be passed to the references subtask to allow it to
        load its schema from disk. Optional, but must be specified if
        ``refSchema`` is not; if both are specified, ``refSchema`` takes
        precedence.
    refSchema : `lsst.afw.table.Schema`, optional
        The schema of the reference catalog, passed to the constructor of the
        references subtask. Optional, but must be specified if ``butler`` is
        not; if both are specified, ``refSchema`` takes precedence.
    **kwds
        Keyword arguments are passed to the supertask constructor.
    """

    ConfigClass = ForcedPhotCoaddConfig
    RunnerClass = ForcedPhotCoaddRunner
    _DefaultName = "forcedPhotCoadd"
    dataPrefix = "deepCoadd_"

    def __init__(self, butler=None, refSchema=None, initInputs=None, **kwds):
        super().__init__(**kwds)

        if initInputs is not None:
            refSchema = initInputs['inputSchema'].schema

        self.makeSubtask("references", butler=butler, schema=refSchema)
        if refSchema is None:
            refSchema = self.references.schema
        self.makeSubtask("measurement", refSchema=refSchema)
        # It is necessary to get the schema internal to the forced measurement task until such a time
        # that the schema is not owned by the measurement task, but is passed in by an external caller
        if self.config.doApCorr:
            self.makeSubtask("applyApCorr", schema=self.measurement.schema)
        self.makeSubtask('catalogCalculation', schema=self.measurement.schema)
        self.outputSchema = lsst.afw.table.SourceCatalog(self.measurement.schema)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        refCatInBand = inputs.pop('refCatInBand')
        inputs['measCat'], inputs['exposureId'] = self.generateMeasCat(inputRefs.exposure.dataId,
                                                                       inputs['exposure'],
                                                                       inputs['refCat'],
                                                                       refCatInBand,
                                                                       inputs['refWcs'],
                                                                       "tract_patch")
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def generateMeasCat(self, exposureDataId, exposure, refCat, refCatInBand, refWcs, idPackerName):
        """Generate a measurement catalog for Gen3.

        Parameters
        ----------
        exposureDataId : `DataId`
            Butler dataId for this exposure.
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure to generate the catalog for.
        refCat : `lsst.afw.table.SourceCatalog`
            Catalog of shapes and positions at which to force photometry.
        refCatInBand : `lsst.afw.table.SourceCatalog`
            Catalog of shapes and position in the band forced photometry is
            currently being performed
        refWcs : `lsst.afw.image.SkyWcs`
            Reference world coordinate system.
        idPackerName : `str`
            Type of ID packer to construct from the registry.

        Returns
        -------
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog of forced sources to measure.
        expId : `int`
            Unique binary id associated with the input exposure

        Raises
        ------
        LookupError
            Raised if a footprint with a given source id was in the reference
            catalog but not in the reference catalog in band (meaning there
            was some sort of mismatch in the two input catalogs)
        """
        expId, expBits = exposureDataId.pack(idPackerName, returnMaxBits=True)
        idFactory = lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

        measCat = self.measurement.generateMeasCat(exposure, refCat, refWcs,
                                                   idFactory=idFactory)
        # attach footprints here, as the attachFootprints method is geared for gen2
        # and is not worth modifying, as this can naturally live inside this method
        for srcRecord in measCat:
            fpRecord = refCatInBand.find(srcRecord.getId())
            if fpRecord is None:
                raise LookupError("Cannot find Footprint for source {}; please check that {} "
                                  "IDs are compatible with reference source IDs"
                                  .format(srcRecord.getId(), self.config.connections.refCatInBand))
            srcRecord.setFootprint(fpRecord.getFootprint())
        return measCat, expId

    def runDataRef(self, dataRef, psfCache=None):
        """Perform forced measurement on a single exposure.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Passed to the ``references`` subtask to obtain the reference WCS,
            the ``getExposure`` method (implemented by derived classes) to
            read the measurment image, and the ``fetchReferences`` method to
            get the exposure and load the reference catalog (see
            :lsst-task`lsst.meas.base.references.CoaddSrcReferencesTask`).
            Refer to derived class documentation for details of the datasets
            and data ID keys which are used.
        psfCache : `int`, optional
            Size of PSF cache, or `None`. The size of the PSF cache can have
            a significant effect upon the runtime for complicated PSF models.

        Notes
        -----
        Sources are generated with ``generateMeasCat`` in the ``measurement``
        subtask. These are passed to ``measurement``'s ``run`` method, which
        fills the source catalog with the forced measurement results. The
        sources are then passed to the ``writeOutputs`` method (implemented by
        derived classes) which writes the outputs.
        """
        refWcs = self.references.getWcs(dataRef)
        exposure = self.getExposure(dataRef)
        if psfCache is not None:
            exposure.getPsf().setCacheCapacity(psfCache)
        refCat = self.fetchReferences(dataRef, exposure)

        measCat = self.measurement.generateMeasCat(exposure, refCat, refWcs,
                                                   idFactory=self.makeIdFactory(dataRef))
        self.log.info("Performing forced measurement on %s" % (dataRef.dataId,))
        self.attachFootprints(measCat, refCat, exposure, refWcs, dataRef)

        exposureId = self.getExposureId(dataRef)

        forcedPhotResult = self.run(measCat, exposure, refCat, refWcs, exposureId=exposureId)

        self.writeOutput(dataRef, forcedPhotResult.measCat)

    def run(self, measCat, exposure, refCat, refWcs, exposureId=None):
        """Perform forced measurement on a single exposure.

        Parameters
        ----------
        measCat : `lsst.afw.table.SourceCatalog`
            The measurement catalog, based on the sources listed in the
            reference catalog.
        exposure : `lsst.afw.image.Exposure`
            The measurement image upon which to perform forced detection.
        refCat : `lsst.afw.table.SourceCatalog`
            The reference catalog of sources to measure.
        refWcs : `lsst.afw.image.SkyWcs`
            The WCS for the references.
        exposureId : `int`
            Optional unique exposureId used for random seed in measurement
            task.

        Returns
        -------
        result : ~`lsst.pipe.base.Struct`
            Structure with fields:

            ``measCat``
                Catalog of forced measurement results
                (`lsst.afw.table.SourceCatalog`).
        """
        self.measurement.run(measCat, exposure, refCat, refWcs, exposureId=exposureId)
        if self.config.doApCorr:
            self.applyApCorr.run(
                catalog=measCat,
                apCorrMap=exposure.getInfo().getApCorrMap()
            )
        self.catalogCalculation.run(measCat)

        return pipeBase.Struct(measCat=measCat)

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
            return self.measurement.attachTransformedFootprints(sources, refCat, exposure, refWcs)

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

    def getExposure(self, dataRef):
        """Read input exposure on which measurement will be performed.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.
        """
        if self.config.hasFakes:
            name = "fakes_" + self.config.coaddName + "Coadd_calexp"
        else:
            name = self.config.coaddName + "Coadd_calexp"

        return dataRef.get(name) if dataRef.datasetExists(name) else None

    def writeOutput(self, dataRef, sources):
        """Write forced source table

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference. The forced_src dataset (with
            self.dataPrefix prepended) is all that will be modified.
        sources : `lsst.afw.table.SourceCatalog`
            Catalog of sources to save.
        """
        dataRef.put(sources, self.dataPrefix + "forced_src", flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)

    def getSchemaCatalogs(self):
        """The schema catalogs that will be used by this task.

        Returns
        -------
        schemaCatalogs : `dict`
            Dictionary mapping dataset type to schema catalog.

        Notes
        -----
        There is only one schema for each type of forced measurement. The
        dataset type for this measurement is defined in the mapper.
        """
        catalog = lsst.afw.table.SourceCatalog(self.measurement.schema)
        catalog.getTable().setMetadata(self.measurement.algMetadata)
        datasetType = self.dataPrefix + "forced_src"
        return {datasetType: catalog}

    def _getConfigName(self):
        # Documented in superclass
        return self.dataPrefix + "forced_config"

    def _getMetadataName(self):
        # Documented in superclass
        return self.dataPrefix + "forced_metadata"

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=lsst.coadd.utils.CoaddDataIdContainer)
        parser.add_argument("--psfCache", type=int, default=100, help="Size of CoaddPsf cache")
        return parser
