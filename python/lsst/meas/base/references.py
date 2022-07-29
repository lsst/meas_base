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

"""
Subtasks for creating the reference catalogs used in forced measurement.
"""

import lsst.geom
import lsst.pex.config
import lsst.pipe.base

__all__ = ("BaseReferencesTask", "CoaddSrcReferencesTask")


class BaseReferencesConfig(lsst.pex.config.Config):
    """Default configuration for reference source selection.
    """

    removePatchOverlaps = lsst.pex.config.Field(
        doc="Only include reference sources for each patch that lie within the patch's inner bbox",
        dtype=bool,
        default=True
    )
    filter = lsst.pex.config.Field(
        doc="Bandpass for reference sources; None indicates chi-squared detections.",
        dtype=str,
        optional=True
    )


class BaseReferencesTask(lsst.pipe.base.Task):
    """Base class for forced photometry subtask that fetches reference sources.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        The schema of the reference catalog.
    butler : `None`.
        Compatibility parameter. Should always be `None`.

    Notes
    -----
    Provides the ``subset`` method, which may be of use to derived
    classes when reimplementing ``fetchInBox``.
    """

    ConfigClass = BaseReferencesConfig
    """Configuration class associated with this task (`lsst.pex.config.Config`).
    """

    def __init__(self, butler=None, schema=None, **kwargs):
        if butler is not None:
            raise ValueError("Gen2 butler parameter is no longer supported.")
        if schema is None:
            raise ValueError("A schema parameter must be provided.")
        lsst.pipe.base.Task.__init__(self, **kwargs)

    def subset(self, sources, bbox, wcs):
        """Filter a list of sources to only those within the bounding box.

        Parameters
        ----------
        sources : iterable of `~lsst.afw.table.SourceRecord`
            Reference sources. May be any Python iterable, including a lazy
            iterator.
        bbox : `lsst.afw.geom.Box2I` or `lsst.afw.geom.Box2D`
            Defines the selection region.
        wcs : `lsst.afw.image.SkyWcs`
            Maps ``bbox`` to sky coordinates.

        Returns
        -------
        sources : iterable of `~lsst.afw.table.SourceRecord`
            Filtered sources. May be any Python iterable, including a lazy
            iterator.

        Notes
        -----
        Instead of filtering sources directly via their positions, we filter
        based on the positions of parent objects, then include or discard all
        children based on their parent's status. This is necessary to support
        replacement with noise in measurement, which requires all child
        sources have their parent present.

        This is not a part of the required `BaseReferencesTask` interface;
        it's a convenience function used in implementing `fetchInBox` that may
        be of use to subclasses.
        """
        boxD = lsst.geom.Box2D(bbox)
        # We're passed an arbitrary iterable, but we need a catalog so we can
        # iterate over parents and then children.
        catalog = lsst.afw.table.SourceCatalog(self.schema)
        catalog.extend(sources)
        # catalog must be sorted by parent ID for lsst.afw.table.getChildren
        # to work
        catalog.sort(lsst.afw.table.SourceTable.getParentKey())
        # Iterate over objects that have no parent.
        parentSources = catalog.getChildren(0)
        skyCoordList = [source.getCoord() for source in parentSources]
        pixelPosList = wcs.skyToPixel(skyCoordList)
        parentList = [parent for parent, pixel in zip(parentSources, pixelPosList) if boxD.contains(pixel)]
        childrenIter = catalog.getChildren((parent.getId() for parent in parentList))
        for parent, children in zip(parentList, childrenIter):
            yield parent
            yield from children


class CoaddSrcReferencesConfig(BaseReferencesTask.ConfigClass):
    """Default configuration for coadd reference source selection.
    """

    coaddName = lsst.pex.config.Field(
        doc="Coadd name: typically one of deep or goodSeeing.",
        dtype=str,
        default="deep",
    )
    skipMissing = lsst.pex.config.Field(
        doc="Silently skip patches where the reference catalog does not exist.",
        dtype=bool,
        default=False
    )

    def validate(self):
        if (self.coaddName == "chiSquared") != (self.filter is None):
            raise lsst.pex.config.FieldValidationError(
                field=CoaddSrcReferencesConfig.coaddName,
                config=self,
                msg="filter may be None if and only if coaddName is chiSquared"
            )


class CoaddSrcReferencesTask(BaseReferencesTask):
    """Select reference sources by loading the “coadd source” dataset directly.

    The name of the dataset to read is generated by appending the
    `datasetSuffix` attribute to the string ``Coadd_``. The dataset is then
    read directly from disk using the Butler.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`, optional
        The schema of the detection catalogs used as input to this one.
    butler : `None`
        Compatibility parameter. Should always be `None`.

    Notes
    -----
    The task will set its own ``self.schema`` attribute to the schema of the
    output merged catalog.
    """

    ConfigClass = CoaddSrcReferencesConfig
    """Configuration class associated with this task (`lsst.pex.config.Config`).
    """

    datasetSuffix = "src"
    """Suffix to append to ``Coadd_`` to generate the dataset name (`str`).
    """

    def __init__(self, butler=None, schema=None, **kwargs):
        BaseReferencesTask.__init__(self, butler=butler, schema=schema, **kwargs)
        self.schema = schema


class MultiBandReferencesConfig(CoaddSrcReferencesTask.ConfigClass):
    """Default configuration for multi-band reference source selection.
    """

    def validate(self):
        if self.filter is not None:
            raise lsst.pex.config.FieldValidationError(
                field=MultiBandReferencesConfig.filter,
                config=self,
                msg="Filter should not be set for the multiband processing scheme")
        # Delegate to ultimate base class, because the direct one has a check we don't want.
        BaseReferencesTask.ConfigClass.validate(self)


class MultiBandReferencesTask(CoaddSrcReferencesTask):
    """Loads references from the multi-band processing scheme.
    """

    ConfigClass = MultiBandReferencesConfig
    datasetSuffix = "ref"  # Documented in superclass
