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

import lsst.pex.config
import lsst.pipe.base
from lsst.utils.logging import PeriodicLogger

import astropy.table

from .baseMeasurement import SimpleBaseMeasurementConfig, SimpleBaseMeasurementTask
from .forcedMeasurement import ForcedPlugin



class SimpleForcedMeasurementConfig(SimpleBaseMeasurementConfig):

    plugins = ForcedPlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_TransformedCentroidFromCoord"
                 "base_PsfFlux",
                 ],
        doc="Plugins to be run and their configuration"
    )
    refCatIdColumn = lsst.pex.config.Field(
        dtype=str,
        default="diaObjectId",
        doc=(
            "Name of the column that provides the object ID from the refCat connection. "
            "measurement.copyColumns['id'] must be set to this value as well."
            "Ignored if refCatStorageClass='SourceCatalog'."
        )
    )
    refCatRaColumn = lsst.pex.config.Field(
        dtype=str,
        default="ra",
        doc=(
            "Name of the column that provides the right ascension (in floating-point degrees) from the "
            "refCat connection. "
            "Ignored if refCatStorageClass='SourceCatalog'."
        )
    )
    refCatDecColumn = lsst.pex.config.Field(
        dtype=str,
        default="dec",
        doc=(
            "Name of the column that provides the declination (in floating-point degrees) from the "
            "refCat connection. "
            "Ignored if refCatStorageClass='SourceCatalog'."
        )
    )
    psfFootprintScaling = lsst.pex.config.Field(
        dtype=float,
        doc="Scaling factor to apply to the PSF shape when footprintSource='psf' (ignored otherwise).",
        default=3.0,
    )

    def setDefaults(self):
        self.slots.centroid = "base_TransformedCentroidFromCoord"
        self.slots.shape = None
        self.slots.apFlux = None
        self.slots.modelFlux = None
        self.slots.psfFlux = "base_PsfFlux"
        self.slots.gaussianFlux = None
        self.slots.calibFlux = None
        self.doReplaceWithNoise = False


class SimpleForcedMeasurementTask(SimpleBaseMeasurementTask):

    ConfigClass = SimpleForcedMeasurementConfig

    def __init__(self, algMetadata=None, **kwds):
        super().__init__(algMetadata=algMetadata, **kwds)
        self.mapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema(), True)
        self.config.slots.setupSchema(self.mapper.editOutputSchema())
        self.initializePlugins(schemaMapper=self.mapper)
        self.addInvalidPsfFlag(self.mapper.editOutputSchema())
        self.schema = self.mapper.getOutputSchema()
        self.schema.checkUnits(parse_strict=True)

    def run(self, table, exposure, refWcs, beginOrder=None, endOrder=None):
        refCat = self.makeMinimalSourceCatalogFromAstropy(table)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measCat.extend(refCat, mapper=self.mapper)
        measCat.setMetadata(self.algMetadata)
        self.attachPsfShapeFootprints(measCat, exposure, scaling=self.config.psfFootprintScaling)
        self.log.info("Performing forced measurement on %d source%s", len(table),
                      "" if len(table) == 1 else "s")
        # Wrap the task logger into a periodic logger.
        periodicLog = PeriodicLogger(self.log)
        for index, (measRecord, refRecord) in enumerate(*zip(measCat, refCat)):
            if measRecord.getFootprint() is None:
                self.log.warning("Skipping object with ID %s that is off the image.", measRecord.getId())
            self.callMeasure(measRecord, exposure, refRecord, refWcs,
                             beginOrder=beginOrder, endOrder=endOrder)
            # Log a message if it has been a while since the last log.
            periodicLog.log("Forced measurement complete for %d parents (and their children) out of %d",
                            index + 1, len(refCat))
        return lsst.pipe.base.Struct(
            table=self.finishOutputTable(table, measCat),
            catalog=measCat,
        )

    def finishOutputTable(self, table, measCat):
        measTable = measCat.asAstropy()
        del measTable["id"]
        del measTable["coord_ra"]
        del measTable["coord_dec"]
        return astropy.table.hstack([table, measTable], join_type="exact"),

    def makeMinimalSourceCatalogFromAstropy(self, table):
        """Create minimal schema SourceCatalog from an Astropy Table.

        The forced measurement subtask expects this as input.

        Parameters
        ----------
        table : `astropy.table.Table`
            Table with locations and ids.

        Returns
        -------
        outputCatalog : `lsst.afw.table.SourceTable`
            Output catalog with minimal schema.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        outputCatalog = lsst.afw.table.SourceCatalog(schema)
        outputCatalog.reserve(len(table))
        # We should make the columns we grab here configurable.
        for objectId, ra, dec in table.iterrows():
            outputRecord = outputCatalog.addNew()
            outputRecord.setId(objectId)
            outputRecord.setCoord(lsst.geom.SpherePoint(ra, dec, lsst.geom.degrees))
        return outputCatalog

    def attachPsfShapeFootprints(self, sources, exposure, scaling=3):
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
        for record in sources:
            # TODO: vectorize this WCS call
            localPoint = wcs.skyToPixel(record.getCoord())
            localIntPoint = lsst.geom.Point2I(localPoint)
            if not bbox.contains(localIntPoint):
                record.setFootprint(None)
                continue
            ellipse = lsst.afw.geom.ellipses.Ellipse(psf.computeShape(localPoint), localPoint)
            ellipse.getCore().scale(scaling)
            spans = lsst.afw.geom.SpanSet.fromShape(ellipse)
            footprint = lsst.afw.detection.Footprint(spans.clippedTo(bbox), bbox)
            footprint.addPeak(localIntPoint.getX(), localIntPoint.getY(),
                              exposure.image._get(localIntPoint, lsst.afw.image.PARENT))
            record.setFootprint(footprint)
