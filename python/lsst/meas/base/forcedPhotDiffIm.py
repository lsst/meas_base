import collections

from lsst.log import Log
import lsst.afw
import lsst.pipe.base
from .forcedPhotCcd import ForcedPhotCcdTask, ForcedPhotCcdConfig, imageOverlapsTract


class ForcedPhotDiffImConfig(ForcedPhotCcdConfig):
    pass


class ForcedPhotDiffImTask(ForcedPhotCcdTask):

    ConfigClass = ForcedPhotDiffImConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotDiffIm"
    dataPrefix = ""

    def _getConfigName(self):
        # Documented in superclass.
        return self.dataPrefix + "forcedPhotDiffIm_config"

    def _getMetadataName(self):
        # Documented in superclass
        return self.dataPrefix + "forcedPhotDiffIm_metadata"

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_diaSrc",
                               help="data ID with raw CCD keys [+ tract optionally], "
                               "e.g. --id visit=12345 ccd=1,2 [tract=0]",
                               ContainerClass=PerTractDiffImDataIdContainer)
        return parser

    def getExposure(self, dataRef):
        """Read input exposure on which measurement will be performed.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.
        """
        return dataRef.get(self.dataPrefix + "deepDiff_differenceExp", immediate=True)

    def writeOutput(self, dataRef, sources):
        """Write forced source table
        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference. The forced_diaSrc dataset (with
            self.dataPrefix prepended) is all that will be modified.
        sources : `lsst.afw.table.SourceCatalog`
            Catalog of sources to save.
        """
        dataRef.put(sources, self.dataPrefix + "forced_diaSrc", flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)


class PerTractDiffImDataIdContainer(lsst.pipe.base.DataIdContainer):
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
        log = Log.getLogger("meas.base.forcedPhotDiffIm.PerTractDiffImDataIdContainer")
        skymap = None
        visitTract = collections.defaultdict(set)   # Set of tracts for each visit
        visitRefs = collections.defaultdict(list)   # List of data references for each visit
        for dataId in self.idList:
            if "tract" not in dataId:
                # Discover which tracts the data overlaps
                log.info("Reading WCS for components of dataId=%s to determine tracts", dict(dataId))
                if skymap is None:
                    skymap = namespace.butler.get(namespace.config.coaddName + "Coadd_skyMap")

                for ref in namespace.butler.subset("deepDiff_differenceExp", dataId=dataId):
                    if not ref.datasetExists("deepDiff_differenceExp"):
                        continue

                    visit = ref.dataId["visit"]
                    visitRefs[visit].append(ref)

                    md = ref.get("deepDiff_differenceExp_md", immediate=True)
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
