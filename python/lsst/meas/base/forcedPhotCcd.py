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

import lsst.pex.config
import lsst.pex.exceptions
from lsst.log import Log
import lsst.pipe.base
import lsst.geom
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.sphgeom

import lsst.pipe.base.connectionTypes as cT
from .forcedPhotImage import ForcedPhotImageTask, ForcedPhotImageConfig
from lsst.meas.base.recalibrateExposure import RecalibrateExposureTask, RecalibrateExposureConnections

__all__ = ("PerTractCcdDataIdContainer", "ForcedPhotCcdConfig", "ForcedPhotCcdTask", "imageOverlapsTract")


class PerTractCcdDataIdContainer(lsst.pipe.base.DataIdContainer):
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
        log = Log.getLogger("meas.base.forcedPhotCcd.PerTractCcdDataIdContainer")
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
        if (not isinstance(e.message, lsst.pex.exceptions.DomainErrorException) and
                not isinstance(e.message, lsst.pex.exceptions.RuntimeErrorException)):
            raise
        return False

    imagePoly = lsst.sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in imageSkyCorners])
    return tractPoly.intersects(imagePoly)  # "intersects" also covers "contains" or "is contained by"


class ForcedPhotCcdConnections(RecalibrateExposureConnections,
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
    )
    refWcs = cT.Input(
        doc="Reference world coordinate system.",
        name="{inputCoaddName}Coadd.wcs",
        storageClass="Wcs",
        dimensions=["abstract_filter", "skymap", "tract", "patch"],
    )
    measCat = cT.Output(
        doc="Output forced photometry catalog.",
        name="forced_src",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector", "skymap", "tract"],
    )


class ForcedPhotCcdConfig(ForcedPhotImageConfig, pipelineConnections=ForcedPhotCcdConnections):
    recalibrate = lsst.pex.config.ConfigurableField(
        target=RecalibrateExposureTask,
        doc="Recalibrate exposure",
    )


class ForcedPhotCcdTask(ForcedPhotImageTask):
    """A command-line driver for performing forced measurement on CCD images.

    Notes
    -----
    This task is a subclass of
    :lsst-task:`lsst.meas.base.forcedPhotImage.ForcedPhotImageTask` which is
    specifically for doing forced measurement on a single CCD exposure, using
    as a reference catalog the detections which were made on overlapping
    coadds.

    The `run` method (inherited from `ForcedPhotImageTask`) takes a
    `~lsst.daf.persistence.ButlerDataRef` argument that corresponds to a single
    CCD. This should contain the data ID keys that correspond to the
    ``forced_src`` dataset (the output dataset for this task), which are
    typically all those used to specify the ``calexp`` dataset (``visit``,
    ``raft``, ``sensor`` for LSST data) as well as a coadd tract. The tract is
    used to look up the appropriate coadd measurement catalogs to use as
    references (e.g. ``deepCoadd_src``; see
    :lsst-task:`lsst.meas.base.references.CoaddSrcReferencesTask` for more
    information). While the tract must be given as part of the dataRef, the
    patches are determined automatically from the bounding box and WCS of the
    calexp to be measured, and the filter used to fetch references is set via
    the ``filter`` option in the configuration of
    :lsst-task:`lsst.meas.base.references.BaseReferencesTask`).

    In addition to the `run` method, `ForcedPhotCcdTask` overrides several
    methods of `ForcedPhotImageTask` to specialize it for single-CCD
    processing, including `~ForcedPhotImageTask.makeIdFactory`,
    `~ForcedPhotImageTask.fetchReferences`, and
    `~ForcedPhotImageTask.getExposure`. None of these should be called
    directly by the user, though it may be useful to override them further in
    subclasses.
    """

    ConfigClass = ForcedPhotCcdConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCcd"
    dataPrefix = ""

    def __init__(self, *args, **kwargs):
        ForcedPhotImageTask.__init__(self, *args, **kwargs)
        self.makeSubtask("recalibrate")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs['refCat'] = self.filterReferences(inputs['exposure'], inputs['refCat'], inputs['refWcs'])
        inputs['measCat'] = self.generateMeasCat(inputRefs.exposure.dataId,
                                                 inputs['exposure'],
                                                 inputs['refCat'], inputs['refWcs'],
                                                 "visit_detector")
        calibs = self.recalibrate.extractCalibs(inputs, single=True)
        inputs["exposure"] = self.recalibrate.run(inputs["exposure"], **calibs)
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def filterReferences(self, exposure, refCat, refWcs):
        """Filter reference catalog so that all sources are within the
        boundaries of the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure to generate the catalog for.
        refCat : `lsst.afw.table.SourceCatalog`
            Catalog of shapes and positions at which to force photometry.
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
        # not contained within the exposure boundaries.
        sources = type(refCat)(refCat.table)
        for record in refCat:
            if expPolygon.contains(record.getCoord().getVector()):
                sources.append(record)
        refCatIdDict = {ref.getId(): ref.getParent() for ref in sources}

        # Step 3: Cull sources that do not have their parent
        # source in the filtered catalog.  Save two copies of each
        # source.
        refSources = type(refCat)(refCat.table)
        for record in refCat:
            if expPolygon.contains(record.getCoord().getVector()):
                recordId = record.getId()
                topId = recordId
                while (topId > 0):
                    if topId in refCatIdDict:
                        topId = refCatIdDict[topId]
                    else:
                        break
                if topId == 0:
                    refSources.append(record)

        # Step 4: Transform source footprints from the reference
        # coordinates to the exposure coordinates.
        for refRecord in refSources:
            refRecord.setFootprint(refRecord.getFootprint().transform(refWcs,
                                                                      expWcs, expRegion))
        # Step 5: Replace reference catalog with filtered source list.
        return refSources

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
        expBits = dataRef.get("ccdExposureId_bits")
        expId = int(dataRef.get("ccdExposureId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

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
                    self.log.warn("Skipping reference %s (child of %s) with bad Footprint",
                                  record.getId(), record.getParent())
                else:
                    self.log.warn("Skipping reference parent %s with bad Footprint", record.getId())
                    badParents.add(record.getId())
            elif record.getParent() not in badParents:
                references.append(record)
        # catalog must be sorted by parent ID for lsst.afw.table.getChildren to work
        references.sort(lsst.afw.table.SourceTable.getParentKey())
        return references

    def getExposure(self, dataRef):
        """Read input exposure for measurement.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.
        """
        exposure = dataRef.get("calexp")
        return self.recalibrate.runDataRef(dataRef, exposure)

    def _getConfigName(self):
        # Documented in superclass.
        return self.dataPrefix + "forcedPhotCcd_config"

    def _getMetadataName(self):
        # Documented in superclass
        return self.dataPrefix + "forcedPhotCcd_metadata"

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID with raw CCD keys [+ tract optionally], "
                               "e.g. --id visit=12345 ccd=1,2 [tract=0]",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser
