#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import collections

import lsst.pex.logging
import lsst.pex.exceptions
import lsst.afw.table
import lsst.afw.image
import lsst.pipe.base
from lsst.geom import convexHull

class PerTractCcdDataIdContainer(lsst.pipe.base.DataIdContainer):
    """A version of lsst.pipe.base.DataIdContainer that combines raw data IDs with a tract.

    Required because we need to add "tract" to the raw data ID keys (defined as whatever we
    use for 'src') when no tract is provided (so that the user is not required to know
    which tracts are spanned by the raw data ID).

    This IdContainer assumes that a calexp is being measured using the detection information,
    a set of reference catalogs, from the set of coadds which intersect with the calexp.
    It needs the calexp id (e.g. visit, raft, sensor), but is also uses the tract to decide
    what set of coadds to use.  The references from the tract whose patches intersect with
    the calexp are used.
    """
    def _addDataRef(self, namespace, dataId, tract):
        """Construct a dataRef based on dataId, but with an added tract key"""
        forcedDataId = dataId.copy()
        forcedDataId['tract'] = tract
        dataRef = namespace.butler.dataRef(datasetType=self.datasetType, dataId=forcedDataId)
        self.refList.append(dataRef)

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")
        log = lsst.pex.logging.Log.getDefaultLog()
        skymap = None
        visitTract = collections.defaultdict(set)   # Set of tracts for each visit
        visitRefs = collections.defaultdict(list)   # List of data references for each visit
        for dataId in self.idList:
            if "tract" not in dataId:
                # Discover which tracts the data overlaps
                log.info("Reading WCS for components of dataId=%s to determine tracts" % (dict(dataId),))
                if skymap is None:
                    skymap = namespace.butler.get(namespace.config.coaddName + "Coadd_skyMap")

                for ref in namespace.butler.subset("calexp", dataId=dataId):
                    if not ref.datasetExists("calexp"):
                        continue

                    visit = ref.dataId["visit"]
                    visitRefs[visit].append(ref)

                    md = ref.get("calexp_md", immediate=True)
                    wcs = lsst.afw.image.makeWcs(md)
                    box = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(0, 0),
                                              lsst.afw.geom.Point2D(md.get("NAXIS1"), md.get("NAXIS2")))
                    # Going with just the nearest tract.  Since we're throwing all tracts for the visit
                    # together, this shouldn't be a problem unless the tracts are much smaller than a CCD.
                    tract = skymap.findTract(wcs.pixelToSky(box.getCenter()))
                    if overlapsTract(tract, wcs, box):
                        visitTract[visit].add(tract.getId())
            else:
                tract = dataId.pop("tract")
                # making a DataRef for src fills out any missing keys and allows us to iterate
                for ref in namespace.butler.subset("src", dataId=dataId):
                    self._addDataRef(namespace, ref.dataId, tract)

        # Ensure all components of a visit are kept together by putting them all in the same set of tracts
        for visit, tractSet in visitTract.iteritems():
            for ref in visitRefs[visit]:
                for tract in tractSet:
                    self._addDataRef(namespace, ref.dataId, tract)
        if visitTract:
            tractCounter = collections.Counter()
            for tractSet in visitTract.itervalues():
                tractCounter.update(tractSet)
            log.info("Number of visits for each tract: %s" % (dict(tractCounter),))


def overlapsTract(tract, imageWcs, imageBox):
    """Return whether the image (specified by Wcs and bounding box) overlaps the tract

    @param tract: TractInfo specifying a tract
    @param imageWcs: Wcs for image
    @param imageBox: Bounding box for image
    @return bool
    """
    tractWcs = tract.getWcs()
    tractCorners = [tractWcs.pixelToSky(lsst.afw.geom.Point2D(coord)).getVector() for
                    coord in tract.getBBox().getCorners()]
    tractPoly = convexHull(tractCorners)

    try:
        imageCorners = [imageWcs.pixelToSky(lsst.afw.geom.Point2D(pix)) for pix in imageBox.getCorners()]
    except lsst.pex.exceptions.LsstCppException, e:
        # Protecting ourselves from awful Wcs solutions in input images
        if (not isinstance(e.message, lsst.pex.exceptions.DomainErrorException) and
            not isinstance(e.message, lsst.pex.exceptions.RuntimeErrorException)):
            raise
        return False

    imagePoly = convexHull([coord.getVector() for coord in imageCorners])
    if imagePoly is None:
        return False
    return tractPoly.intersects(imagePoly) # "intersects" also covers "contains" or "is contained by"
