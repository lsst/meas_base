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

import collections
import logging
import pandas as pd
import numpy as np

import lsst.pex.config
import lsst.pex.exceptions
import lsst.pipe.base
import lsst.geom
import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.sphgeom

from lsst.obs.base import ExposureIdInfo
from lsst.pipe.base import PipelineTaskConnections
import lsst.pipe.base.connectionTypes as cT

import lsst.pipe.base as pipeBase
from lsst.skymap import BaseSkyMap

from .references import MultiBandReferencesTask
from .forcedMeasurement import ForcedMeasurementTask
from .applyApCorr import ApplyApCorrTask
from .catalogCalculation import CatalogCalculationTask

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("PerTractCcdDataIdContainer", "ForcedPhotCcdConfig", "ForcedPhotCcdTask", "imageOverlapsTract",
           "ForcedPhotCcdFromDataFrameTask", "ForcedPhotCcdFromDataFrameConfig")


class PerTractCcdDataIdContainer(pipeBase.DataIdContainer):
    """A data ID container which combines raw data IDs with a tract.

    Notes
    -----
    Required because we need to add "tract" to the raw data ID keys (defined as
    whatever we use for ``src``) when no tract is provided (so that the user is
    not required to know which tracts are spanned by the raw data ID).

    This subclass of `~lsst.pipe.base.DataIdContainer` assumes that a calexp is
    being measured using the detection information, a set of reference
    catalogs, from the set of coadds which intersect with the calexp.  It needs
    the calexp id (e.g.  visit, raft, sensor), but is also uses the tract to
    decide what set of coadds to use.  The references from the tract whose
    patches intersect with the calexp are used.
    """

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")
        log = logging.getLogger("meas.base.forcedPhotCcd.PerTractCcdDataIdContainer")
        skymap = None
        visitTract = collections.defaultdict(set)   # Set of tracts for each visit
        visitRefs = collections.defaultdict(list)   # List of data references for each visit
        for dataId in self.idList:
            if "tract" not in dataId:
                # Discover which tracts the data overlaps
                log.info("Reading WCS for components of dataId=%s to determine tracts", dict(dataId))
                if skymap is None:
                    skymap = namespace.butler.get(namespace.config.coaddName + "Coadd_skyMap")

                for ref in namespace.butler.subset("calexp", dataId=dataId):
                    if not ref.datasetExists("calexp"):
                        continue

                    visit = ref.dataId["visit"]
                    visitRefs[visit].append(ref)

                    md = ref.get("calexp_md", immediate=True)
                    wcs = lsst.afw.geom.makeSkyWcs(md)
                    box = lsst.geom.Box2D(lsst.afw.image.bboxFromMetadata(md))
                    # Going with just the nearest tract.  Since we're throwing all tracts for the visit
                    # together, this shouldn't be a problem unless the tracts are much smaller than a CCD.
                    tract = skymap.findTract(wcs.pixelToSky(box.getCenter()))
                    if imageOverlapsTract(tract, wcs, box):
                        visitTract[visit].add(tract.getId())
            else:
                self.refList.extend(ref for ref in namespace.butler.subset(self.datasetType, dataId=dataId))

        # Ensure all components of a visit are kept together by putting them all in the same set of tracts
        for visit, tractSet in visitTract.items():
            for ref in visitRefs[visit]:
                for tract in tractSet:
                    self.refList.append(namespace.butler.dataRef(datasetType=self.datasetType,
                                                                 dataId=ref.dataId, tract=tract))
        if visitTract:
            tractCounter = collections.Counter()
            for tractSet in visitTract.values():
                tractCounter.update(tractSet)
            log.info("Number of visits for each tract: %s", dict(tractCounter))


def imageOverlapsTract(tract, imageWcs, imageBox):
    """Return whether the given bounding box overlaps the tract given a WCS.

    Parameters
    ----------
    tract : `lsst.skymap.TractInfo`
        TractInfo specifying a tract.
    imageWcs : `lsst.afw.geom.SkyWcs`
        World coordinate system for the image.
    imageBox : `lsst.geom.Box2I`
        Bounding box for the image.

    Returns
    -------
    overlap : `bool`
        `True` if the bounding box overlaps the tract; `False` otherwise.
    """
    tractPoly = tract.getOuterSkyPolygon()

    imagePixelCorners = lsst.geom.Box2D(imageBox).getCorners()
    try:
        imageSkyCorners = imageWcs.pixelToSky(imagePixelCorners)
    except lsst.pex.exceptions.LsstCppException as e:
        # Protecting ourselves from awful Wcs solutions in input images
        if (not isinstance(e.message, lsst.pex.exceptions.DomainErrorException)
                and not isinstance(e.message, lsst.pex.exceptions.RuntimeErrorException)):
            raise
        return False

    imagePoly = lsst.sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in imageSkyCorners])
    return tractPoly.intersects(imagePoly)  # "intersects" also covers "contains" or "is contained by"


class ForcedPhotCcdConnections(PipelineTaskConnections,
                               dimensions=("instrument", "visit", "detector", "skymap", "tract"),
                               defaultTemplates={"inputCoaddName": "deep",
                                                 "inputName": "calexp"}):
    inputSchema = cT.InitInput(
        doc="Schema for the input measurement catalogs.",
        name="{inputCoaddName}Coadd_ref_schema",
        storageClass="SourceCatalog",
    )
    outputSchema = cT.InitOutput(
        doc="Schema for the output forced measurement catalogs.",
        name="forced_src_schema",
        storageClass="SourceCatalog",
    )
    exposure = cT.Input(
        doc="Input exposure to perform photometry on.",
        name="{inputName}",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )
    refCat = cT.Input(
        doc="Catalog of shapes and positions at which to force photometry.",
        name="{inputCoaddName}Coadd_ref",
        storageClass="SourceCatalog",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
    )
    skyMap = cT.Input(
        doc="SkyMap dataset that defines the coordinate system of the reference catalog.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=["skymap"],
    )
    measCat = cT.Output(
        doc="Output forced photometry catalog.",
        name="forced_src",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector", "skymap", "tract"],
    )


class ForcedPhotCcdConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=ForcedPhotCcdConnections):
    """Config class for forced measurement driver task."""
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
    doApplyUberCal = lsst.pex.config.Field(
        dtype=bool,
        doc="Apply meas_mosaic ubercal results to input calexps?",
        default=False,
        deprecated="Deprecated by DM-23352; use doApplyExternalPhotoCalib and doApplyExternalSkyWcs instead",
    )
    doApplyExternalPhotoCalib = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("Whether to apply external photometric calibration via an "
             "`lsst.afw.image.PhotoCalib` object. Uses the "
             "``externalPhotoCalibName`` field to determine which calibration "
             "to load."),
    )
    doApplyExternalSkyWcs = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("Whether to apply external astrometric calibration via an "
             "`lsst.afw.geom.SkyWcs` object. Uses ``externalSkyWcsName`` "
             "field to determine which calibration to load."),
    )
    doApplySkyCorr = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Apply sky correction?",
    )
    includePhotoCalibVar = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Add photometric calibration variance to warp variance plane?",
    )
    externalPhotoCalibName = lsst.pex.config.ChoiceField(
        dtype=str,
        doc=("Type of external PhotoCalib if ``doApplyExternalPhotoCalib`` is True. "
             "Unused for Gen3 middleware."),
        default="jointcal",
        allowed={
            "jointcal": "Use jointcal_photoCalib",
            "fgcm": "Use fgcm_photoCalib",
            "fgcm_tract": "Use fgcm_tract_photoCalib"
        },
    )
    externalSkyWcsName = lsst.pex.config.ChoiceField(
        dtype=str,
        doc="Type of external SkyWcs if ``doApplyExternalSkyWcs`` is True. Unused for Gen3 middleware.",
        default="jointcal",
        allowed={
            "jointcal": "Use jointcal_wcs"
        },
    )
    footprintSource = lsst.pex.config.ChoiceField(
        dtype=str,
        doc="Where to obtain footprints to install in the measurement catalog, prior to measurement.",
        allowed={
            "transformed": "Transform footprints from the reference catalog (downgrades HeavyFootprints).",
            "psf": ("Use the scaled shape of the PSF at the position of each source (does not generate "
                    "HeavyFootprints)."),
        },
        optional=True,
        default="transformed",
    )
    psfFootprintScaling = lsst.pex.config.Field(
        dtype=float,
        doc="Scaling factor to apply to the PSF shape when footprintSource='psf' (ignored otherwise).",
        default=3.0,
    )

    def setDefaults(self):
        # Docstring inherited.
        # Make catalogCalculation a no-op by default as no modelFlux is setup by default in
        # ForcedMeasurementTask
        super().setDefaults()
        self.measurement.plugins.names |= ['base_LocalPhotoCalib', 'base_LocalWcs']
        self.catalogCalculation.plugins.names = []


class ForcedPhotCcdTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    """A command-line driver for performing forced measurement on CCD images.

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

    Notes
    -----
    The `runDataRef` method takes a `~lsst.daf.persistence.ButlerDataRef` argument
    that corresponds to a single CCD. This should contain the data ID keys that
    correspond to the ``forced_src`` dataset (the output dataset for this
    task), which are typically all those used to specify the ``calexp`` dataset
    (``visit``, ``raft``, ``sensor`` for LSST data) as well as a coadd tract.
    The tract is used to look up the appropriate coadd measurement catalogs to
    use as references (e.g. ``deepCoadd_src``; see
    :lsst-task:`lsst.meas.base.references.CoaddSrcReferencesTask` for more
    information). While the tract must be given as part of the dataRef, the
    patches are determined automatically from the bounding box and WCS of the
    calexp to be measured, and the filter used to fetch references is set via
    the ``filter`` option in the configuration of
    :lsst-task:`lsst.meas.base.references.BaseReferencesTask`).
    """

    ConfigClass = ForcedPhotCcdConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCcd"
    dataPrefix = ""

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

        tract = butlerQC.quantum.dataId['tract']
        skyMap = inputs.pop("skyMap")
        inputs['refWcs'] = skyMap[tract].getWcs()

        inputs['refCat'] = self.mergeAndFilterReferences(inputs['exposure'], inputs['refCat'],
                                                         inputs['refWcs'])

        inputs['measCat'], inputs['exposureId'] = self.generateMeasCat(inputRefs.exposure.dataId,
                                                                       inputs['exposure'],
                                                                       inputs['refCat'], inputs['refWcs'],
                                                                       "visit_detector")
        self.attachFootprints(inputs['measCat'], inputs['refCat'], inputs['exposure'], inputs['refWcs'])
        # TODO: apply external calibrations (DM-17062)
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def mergeAndFilterReferences(self, exposure, refCats, refWcs):
        """Filter reference catalog so that all sources are within the
        boundaries of the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure to generate the catalog for.
        refCats : sequence of `lsst.daf.butler.DeferredDatasetHandle`
            Handles for catalogs of shapes and positions at which to force
            photometry.
        refWcs : `lsst.afw.image.SkyWcs`
            Reference world coordinate system.

        Returns
        -------
        refSources : `lsst.afw.table.SourceCatalog`
            Filtered catalog of forced sources to measure.

        Notes
        -----
        Filtering the reference catalog is currently handled by Gen2
        specific methods.  To function for Gen3, this method copies
        code segments to do the filtering and transformation.  The
        majority of this code is based on the methods of
        lsst.meas.algorithms.loadReferenceObjects.ReferenceObjectLoader

        """

        # Step 1: Determine bounds of the exposure photometry will
        # be performed on.
        expWcs = exposure.getWcs()
        expRegion = exposure.getBBox(lsst.afw.image.PARENT)
        expBBox = lsst.geom.Box2D(expRegion)
        expBoxCorners = expBBox.getCorners()
        expSkyCorners = [expWcs.pixelToSky(corner).getVector() for
                         corner in expBoxCorners]
        expPolygon = lsst.sphgeom.ConvexPolygon(expSkyCorners)

        # Step 2: Filter out reference catalog sources that are
        # not contained within the exposure boundaries, or whose
        # parents are not within the exposure boundaries.  Note
        # that within a single input refCat, the parents always
        # appear before the children.
        mergedRefCat = None
        for refCat in refCats:
            refCat = refCat.get()
            if mergedRefCat is None:
                mergedRefCat = lsst.afw.table.SourceCatalog(refCat.table)
            containedIds = {0}  # zero as a parent ID means "this is a parent"
            for record in refCat:
                if expPolygon.contains(record.getCoord().getVector()) and record.getParent() in containedIds:
                    record.setFootprint(record.getFootprint())
                    mergedRefCat.append(record)
                    containedIds.add(record.getId())
        if mergedRefCat is None:
            raise RuntimeError("No reference objects for forced photometry.")
        mergedRefCat.sort(lsst.afw.table.SourceTable.getParentKey())
        return mergedRefCat

    def generateMeasCat(self, exposureDataId, exposure, refCat, refWcs, idPackerName):
        """Generate a measurement catalog for Gen3.

        Parameters
        ----------
        exposureDataId : `DataId`
            Butler dataId for this exposure.
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure to generate the catalog for.
        refCat : `lsst.afw.table.SourceCatalog`
            Catalog of shapes and positions at which to force photometry.
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
        """
        exposureIdInfo = ExposureIdInfo.fromDataId(exposureDataId, idPackerName)
        idFactory = exposureIdInfo.makeSourceIdFactory()

        measCat = self.measurement.generateMeasCat(exposure, refCat, refWcs,
                                                   idFactory=idFactory)
        return measCat, exposureIdInfo.expId

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
            exposure.getPsf().setCacheSize(psfCache)
        refCat = self.fetchReferences(dataRef, exposure)
        measCat = self.measurement.generateMeasCat(exposure, refCat, refWcs,
                                                   idFactory=self.makeIdFactory(dataRef))
        self.log.info("Performing forced measurement on %s", dataRef.dataId)
        self.attachFootprints(measCat, refCat, exposure, refWcs)

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
        result : `lsst.pipe.base.Struct`
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
            Butler data reference. The ``ccdExposureId_bits`` and
            ``ccdExposureId`` datasets are accessed. The data ID must have the
            keys that correspond to ``ccdExposureId``, which are generally the
            same as those that correspond to ``calexp`` (``visit``, ``raft``,
            ``sensor`` for LSST data).
        """
        exposureIdInfo = ExposureIdInfo(int(dataRef.get("ccdExposureId")), dataRef.get("ccdExposureId_bits"))
        return exposureIdInfo.makeSourceIdFactory()

    def getExposureId(self, dataRef):
        return int(dataRef.get("ccdExposureId", immediate=True))

    def fetchReferences(self, dataRef, exposure):
        """Get sources that overlap the exposure.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference corresponding to the image to be measured;
            should have ``tract``, ``patch``, and ``filter`` keys.
        exposure : `lsst.afw.image.Exposure`
            The image to be measured (used only to obtain a WCS and bounding
            box).

        Returns
        -------
        referencs : `lsst.afw.table.SourceCatalog`
            Catalog of sources that overlap the exposure

        Notes
        -----
        The returned catalog is sorted by ID and guarantees that all included
        children have their parent included and that all Footprints are valid.

        All work is delegated to the references subtask; see
        :lsst-task:`lsst.meas.base.references.CoaddSrcReferencesTask`
        for information about the default behavior.
        """
        references = lsst.afw.table.SourceCatalog(self.references.schema)
        badParents = set()
        unfiltered = self.references.fetchInBox(dataRef, exposure.getBBox(), exposure.getWcs())
        for record in unfiltered:
            if record.getFootprint() is None or record.getFootprint().getArea() == 0:
                if record.getParent() != 0:
                    self.log.warning("Skipping reference %s (child of %s) with bad Footprint",
                                     record.getId(), record.getParent())
                else:
                    self.log.warning("Skipping reference parent %s with bad Footprint", record.getId())
                    badParents.add(record.getId())
            elif record.getParent() not in badParents:
                references.append(record)
        # catalog must be sorted by parent ID for lsst.afw.table.getChildren to work
        references.sort(lsst.afw.table.SourceTable.getParentKey())
        return references

    def attachFootprints(self, sources, refCat, exposure, refWcs):
        r"""Attach footprints to blank sources prior to measurements.

        Notes
        -----
        `~lsst.afw.detection.Footprint`\ s for forced photometry must be in the
        pixel coordinate system of the image being measured, while the actual
        detections may start out in a different coordinate system.

        Subclasses of this class may implement this method to define how
        those `~lsst.afw.detection.Footprint`\ s should be generated.

        This default implementation transforms depends on the
        ``footprintSource`` configuration parameter.
        """
        if self.config.footprintSource == "transformed":
            return self.measurement.attachTransformedFootprints(sources, refCat, exposure, refWcs)
        elif self.config.footprintSource == "psf":
            return self.measurement.attachPsfShapeFootprints(sources, exposure,
                                                             scaling=self.config.psfFootprintScaling)

    def getExposure(self, dataRef):
        """Read input exposure for measurement.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.
        """
        exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True)

        if self.config.doApplyExternalPhotoCalib:
            source = f"{self.config.externalPhotoCalibName}_photoCalib"
            self.log.info("Applying external photoCalib from %s", source)
            photoCalib = dataRef.get(source)
            exposure.setPhotoCalib(photoCalib)  # No need for calibrateImage; having the photoCalib suffices

        if self.config.doApplyExternalSkyWcs:
            source = f"{self.config.externalSkyWcsName}_wcs"
            self.log.info("Applying external skyWcs from %s", source)
            skyWcs = dataRef.get(source)
            exposure.setWcs(skyWcs)

        if self.config.doApplySkyCorr:
            self.log.info("Apply sky correction")
            skyCorr = dataRef.get("skyCorr")
            exposure.maskedImage -= skyCorr.getImage()

        return exposure

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
        # Documented in superclass.
        return self.dataPrefix + "forcedPhotCcd_config"

    def _getMetadataName(self):
        # Documented in superclass
        return self.dataPrefix + "forcedPhotCcd_metadata"

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID with raw CCD keys [+ tract optionally], "
                               "e.g. --id visit=12345 ccd=1,2 [tract=0]",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser


class ForcedPhotCcdFromDataFrameConnections(PipelineTaskConnections,
                                            dimensions=("instrument", "visit", "detector", "skymap", "tract"),
                                            defaultTemplates={"inputCoaddName": "goodSeeing",
                                                              "inputName": "calexp"}):
    refCat = cT.Input(
        doc="Catalog of positions at which to force photometry.",
        name="{inputCoaddName}Diff_fullDiaObjTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
    )
    exposure = cT.Input(
        doc="Input exposure to perform photometry on.",
        name="{inputName}",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )
    measCat = cT.Output(
        doc="Output forced photometry catalog.",
        name="forced_src_diaObject",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector", "skymap", "tract"],
    )
    outputSchema = cT.InitOutput(
        doc="Schema for the output forced measurement catalogs.",
        name="forced_src_diaObject_schema",
        storageClass="SourceCatalog",
    )


class ForcedPhotCcdFromDataFrameConfig(ForcedPhotCcdConfig,
                                       pipelineConnections=ForcedPhotCcdFromDataFrameConnections):
    def setDefaults(self):
        super().setDefaults()
        self.footprintSource = "psf"
        self.measurement.doReplaceWithNoise = False
        self.measurement.plugins = ["base_TransformedCentroidFromCoord", "base_PsfFlux", "base_PixelFlags"]
        self.measurement.copyColumns = {'id': 'diaObjectId', 'coord_ra': 'coord_ra', 'coord_dec': 'coord_dec'}
        self.measurement.slots.centroid = "base_TransformedCentroidFromCoord"
        self.measurement.slots.psfFlux = "base_PsfFlux"
        self.measurement.slots.shape = None

    def validate(self):
        super().validate()
        if self.footprintSource == "transformed":
            raise ValueError("Cannot transform footprints from reference catalog, "
                             "because DataFrames can't hold footprints.")


class ForcedPhotCcdFromDataFrameTask(ForcedPhotCcdTask):
    """Force Photometry on a per-detector exposure with coords from a DataFrame

    Uses input from a DataFrame instead of SourceCatalog
    like the base class ForcedPhotCcd does.
    Writes out a SourceCatalog so that the downstream
    WriteForcedSourceTableTask can be reused with output from this Task.
    """
    _DefaultName = "forcedPhotCcdFromDataFrame"
    ConfigClass = ForcedPhotCcdFromDataFrameConfig

    def __init__(self, butler=None, refSchema=None, initInputs=None, **kwds):
        # Parent's init assumes that we have a reference schema; Cannot reuse
        pipeBase.PipelineTask.__init__(self, **kwds)

        self.makeSubtask("measurement", refSchema=lsst.afw.table.SourceTable.makeMinimalSchema())

        if self.config.doApCorr:
            self.makeSubtask("applyApCorr", schema=self.measurement.schema)
        self.makeSubtask('catalogCalculation', schema=self.measurement.schema)
        self.outputSchema = lsst.afw.table.SourceCatalog(self.measurement.schema)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        self.log.info("Filtering ref cats: %s", ','.join([str(i.dataId) for i in inputs['refCat']]))
        refCat = self.df2RefCat([i.get(parameters={"columns": ['diaObjectId', 'ra', 'decl']})
                                 for i in inputs['refCat']],
                                inputs['exposure'].getBBox(), inputs['exposure'].getWcs())
        inputs['refCat'] = refCat
        inputs['refWcs'] = inputs['exposure'].getWcs()
        inputs['measCat'], inputs['exposureId'] = self.generateMeasCat(inputRefs.exposure.dataId,
                                                                       inputs['exposure'], inputs['refCat'],
                                                                       inputs['refWcs'],
                                                                       "visit_detector")
        self.attachFootprints(inputs["measCat"], inputs["refCat"], inputs["exposure"], inputs["refWcs"])
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def df2RefCat(self, dfList, exposureBBox, exposureWcs):
        """Convert list of DataFrames to reference catalog

        Concatenate list of DataFrames presumably from multiple patches and
        downselect rows that overlap the exposureBBox using the exposureWcs.

        Parameters
        ----------
        dfList : `list` of `pandas.DataFrame`
            Each element containst diaObjects with ra/decl position in degrees
            Columns 'diaObjectId', 'ra', 'decl' are expected
        exposureBBox :   `lsst.geom.Box2I`
            Bounding box on which to select rows that overlap
        exposureWcs : `lsst.afw.geom.SkyWcs`
            World coordinate system to convert sky coords in ref cat to
            pixel coords with which to compare with exposureBBox

        Returns
        -------
        refCat : `lsst.afw.table.SourceTable`
            Source Catalog with minimal schema that overlaps exposureBBox
        """
        df = pd.concat(dfList)
        # translate ra/decl coords in dataframe to detector pixel coords
        # to down select rows that overlap the detector bbox
        mapping = exposureWcs.getTransform().getMapping()
        x, y = mapping.applyInverse(np.array(df[['ra', 'decl']].values*2*np.pi/360).T)
        inBBox = exposureBBox.contains(x, y)
        refCat = self.df2SourceCat(df[inBBox])
        return refCat

    def df2SourceCat(self, df):
        """Create minimal schema SourceCatalog from a pandas DataFrame.

        The forced measurement subtask expects this as input.

        Parameters
        ----------
        df : `pandas.DataFrame`
            DiaObjects with locations and ids.

        Returns
        -------
        outputCatalog : `lsst.afw.table.SourceTable`
            Output catalog with minimal schema.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        outputCatalog = lsst.afw.table.SourceCatalog(schema)
        outputCatalog.reserve(len(df))

        for id, diaObjectId, ra, decl in df[['diaObjectId', 'ra', 'decl']].itertuples():
            outputRecord = outputCatalog.addNew()
            outputRecord.setId(diaObjectId)
            outputRecord.setCoord(lsst.geom.SpherePoint(ra, decl, lsst.geom.degrees))
        return outputCatalog
