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

__all__ = ("SimpleForcedMeasurementConfig", "SimpleForcedMeasurementTask")

import numpy as np

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.daf.base
import lsst.geom
import lsst.pex.config
import lsst.pipe.base
from lsst.utils.logging import PeriodicLogger

from .baseMeasurement import SimpleBaseMeasurementConfig, SimpleBaseMeasurementTask
from .forcedMeasurement import ForcedPlugin


class SimpleForcedMeasurementConfig(SimpleBaseMeasurementConfig):
    """Config class for SimpleForcedMeasurementTask."""
    plugins = ForcedPlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_TransformedCentroidFromCoord",
                 "base_PsfFlux",
                 ],
        doc="Plugins to be run and their configuration"
    )
    psfFootprintScaling = lsst.pex.config.Field(
        dtype=float,
        doc="Scaling factor to apply to the PSF shape when footprintSource='psf' (ignored otherwise).",
        default=3.0,
    )
    copyColumns = lsst.pex.config.DictField(
        keytype=str, itemtype=str, doc="Mapping of reference columns to source columns",
        default={"id": "objectId", "parent": "parentObjectId",
                 "coord_ra": "coord_ra", "coord_dec": "coord_dec"}
    )
    checkUnitsParseStrict = lsst.pex.config.Field(
        doc="Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'",
        dtype=str,
        default="raise",
    )

    def setDefaults(self):
        self.slots.centroid = "base_TransformedCentroidFromCoord"
        self.slots.shape = None
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = "base_PsfFlux"
        self.slots.gaussianFlux = None
        self.slots.calibFlux = None


class SimpleForcedMeasurementTask(SimpleBaseMeasurementTask):
    """Measure sources on an image using a simple forced measurement algorithm.

    This differes from ForcedMeasurmentTask in that it uses a PSF-based
    footprint for every source so it does not need to transform footprints.

    Parameters
    ----------
    algMetadata : `lsst.daf.base.PropertyList` or `None`
        Will be updated in place to to record information about each
        algorithm. An empty `~lsst.daf.base.PropertyList` will be created if
        `None`.
    **kwds
        Keyword arguments are passed to the supertask constructor.
    """
    ConfigClass = SimpleForcedMeasurementConfig

    def __init__(self, refSchema, algMetadata: lsst.daf.base.PropertyList = None, **kwds):
        super().__init__(algMetadata=algMetadata, **kwds)
        self.mapper = lsst.afw.table.SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema(), False)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        for refName, targetName in self.config.copyColumns.items():
            refItem = refSchema.find(refName)
            self.mapper.addMapping(refItem.key, targetName)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        self.initializePlugins(schemaMapper=self.mapper)
        self.addInvalidPsfFlag(self.mapper.editOutputSchema())
        self.schema = self.mapper.getOutputSchema()
        self.schema.checkUnits(parse_strict=self.config.checkUnitsParseStrict)

    def run(
        self,
        refCat: lsst.afw.table.SourceCatalog,
        measCat: lsst.afw.table.SourceCatalog,
        exposure: lsst.afw.image.Exposure,
        refWcs: lsst.afw.geom.SkyWcs,
        beginOrder: int | None = None,
        endOrder: int | None = None,
    ) -> None:
        """Perform forced measurement.

        Parameters
        ----------
        refCat : `lsst.afw.table.SourceCatalog`
            Catalog with locations and ids of sources to measure.
        measCat : `lsst.afw.table.SourceCatalog`
            Catalog that measurements are made on.
        exposure : `lsst.afw.image.exposureF`
            Image to be measured. Must have at least a `lsst.afw.geom.SkyWcs`
            attached.
        refWcs : `lsst.afw.geom.SkyWcs`
            Defines the X,Y coordinate system of ``refCat``.
        beginOrder : `int`, optional
            Beginning execution order (inclusive). Algorithms with
            ``executionOrder`` < ``beginOrder`` are not executed. `None` for no limit.
        endOrder : `int`, optional
            Ending execution order (exclusive). Algorithms with
            ``executionOrder`` >= ``endOrder`` are not executed. `None` for no limit.
        idFactory : `lsst.afw.table.IdFactory`, optional
            Factory for creating IDs for sources.
        """
        self._attachPsfShapeFootprints(measCat, exposure, scaling=self.config.psfFootprintScaling)
        self.log.info("Performing forced measurement on %d source%s", len(refCat),
                      "" if len(refCat) == 1 else "s")
        # Wrap the task logger into a periodic logger.
        periodicLog = PeriodicLogger(self.log)

        for index in range(len(refCat)):
            measRecord = measCat[index]
            refRecord = refCat[index]
            if measRecord.getFootprint() is None:
                self.log.warning("Skipping object with ID %s that is off the image.", measRecord.getId())
            self.callMeasure(measRecord, exposure, refRecord, refWcs,
                             beginOrder=beginOrder, endOrder=endOrder)
            # Log a message if it has been a while since the last log.
            periodicLog.log("Forced measurement complete for %d parents (and their children) out of %d",
                            index + 1, len(refCat))

    def _attachPsfShapeFootprints(self, sources, exposure, scaling=3):
        """Attach Footprints to blank sources prior to measurement, by
        creating elliptical Footprints from the PSF moments.

        Parameters
        ----------
        sources : `lsst.afw.table.SourceCatalog`
            Blank catalog (with all rows and columns, but values other than
            ``coord_ra``, ``coord_dec`` unpopulated).
            to which footprints should be attached.
        exposure : `lsst.afw.image.Exposure`
            Image object from which peak values and the PSF are obtained.
        scaling : `int`, optional
            Scaling factor to apply to the PSF second-moments ellipse in order
            to determine the footprint boundary.
        """
        psf = exposure.getPsf()
        if psf is None:
            raise RuntimeError("Cannot construct Footprints from PSF shape without a PSF.")
        bbox = exposure.getBBox()
        wcs = exposure.getWcs()
        # This will always be coord_ra, coord_dec since we converted the
        # astropy table into a schema and schema is always coord_ra, coord_dec.
        x, y = wcs.skyToPixelArray(sources["coord_ra"], sources["coord_dec"], degrees=False)
        inBBox = np.atleast_1d(lsst.geom.Box2D(bbox).contains(x, y))
        for idx, record in enumerate(sources):
            localPoint = lsst.geom.Point2D(x[idx], y[idx])
            localIntPoint = lsst.geom.Point2I(localPoint)
            if not inBBox[idx]:
                record.setFootprint(None)
                continue
            ellipse = lsst.afw.geom.ellipses.Ellipse(psf.computeShape(localPoint), localPoint)
            ellipse.getCore().scale(scaling)
            spans = lsst.afw.geom.SpanSet.fromShape(ellipse)
            footprint = lsst.afw.detection.Footprint(spans.clippedTo(bbox), bbox)
            footprint.addPeak(localIntPoint.getX(), localIntPoint.getY(),
                              exposure.image._get(localIntPoint, lsst.afw.image.PARENT))
            record.setFootprint(footprint)
