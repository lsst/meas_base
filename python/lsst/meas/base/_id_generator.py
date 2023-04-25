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

from __future__ import annotations

__all__ = (
    "IdGenerator",
    "FullIdGenerator",
    "BaseIdGeneratorConfig",
    "DetectorExposureIdGeneratorConfig",
    "DetectorVisitIdGeneratorConfig",
    "SkyMapIdGeneratorConfig",
)

import dataclasses
from typing import Any, Callable

import numpy as np
from lsst.afw.table import IdFactory, Schema, SourceCatalog, SourceTable
from lsst.daf.butler import DataCoordinate, DimensionPacker
from lsst.obs.base import ExposureIdInfo
from lsst.pex.config import Config, ConfigField, Field
from lsst.pipe.base import Instrument
from lsst.skymap.packers import SkyMapDimensionPacker

DEFAULT_RELEASE_ID = 0
"""Default release ID to embed in catalog IDs.

This can be changed globally to avoid having to override individual task
configs to set the release ID.
"""

DEFAULT_N_RELEASES = 1  # 1 means don't reserve space for releases.
"""Default number of releases to reserve space for in catalog IDs."""


class BaseIdGeneratorConfig(Config):
    """Base class for configuration of `IdGenerator` instances.

    This class is abstract (it cannot use `abc.ABCMeta` due to a metaclass
    conflict), and it should mostly be considered an implementation detail
    of how the attributes it defines are included in its concrete derived
    classes.  Derived classes must implemented `_make_dimension_packer`.

    See `IdGenerator` for usage.
    """

    release_id = Field(
        doc=(
            "Identifier for a data release or other version to embed in generated IDs. "
            "Zero is reserved for IDs with no embedded release identifier."
        ),
        dtype=int,
        default=DEFAULT_RELEASE_ID,
        check=lambda x: x >= 0,
    )

    n_releases = Field(
        doc=(
            "Number of (contiguous, starting from zero) `release_id` values to reserve space for. "
            "One (not zero) is used to reserve no space."
        ),
        dtype=int,
        default=DEFAULT_N_RELEASES,
        check=lambda x: x > 0,
    )

    @classmethod
    def make_field(
        cls, doc="Configuration for how to generate catalog IDs from data IDs."
    ):
        """Return a config field that holds an instance of this class.

        Parameters
        ----------
        doc : `str`, optional
            Documentation for the config field.  As this configuration almost
            always plays the same role in any parent config, the default is
            usually fine.

        Returns
        -------
        field : `lsst.pex.config.ConfigField`
            New config field for instances of this class.

        Notes
        -----
        This method is provided as a convenience to reduce boilerplate
        downstream: it typically saves an import or two, and it allows the same
        usually-appropriate docstring to be reused instead of rewritten each
        time.  It does not need to be used in order to use this config class.
        """
        return ConfigField(doc, dtype=cls)

    def apply(self, data_id: DataCoordinate, **kwargs: Any) -> IdGenerator:
        """Construct an `IdGenerator` instance from this configuration.

        Parameters
        ----------
        data_id : `DataCoordinate`
            The data ID the `IdGenerator` will embed into all IDs.  This
            generally must be a fully-expanded data ID (i.e. have dimension
            records attached), that identifies the "instrument" or "skymap"
            dimension, though this requirement may be relaxed for certain
            dimension packer types.
        **kwargs
            Additional keyword arguments are interpreted as dimension value
            pairs to include in the data ID.  This may be used to provide
            constraints on dimensions for which records are not available.

        Returns
        -------
        id_generator : `IdGenerator`
            Object that generates integer IDs for catalogs and their rows by
            embedding the given data ID and a configurably-optional release ID.

        Notes
        -----
        This method is called `apply` for consistency with the pattern of using
        `lsst.pex.config.ConfigurableField` and `lsst.pex.config.RegistryField`
        to construct the objects whose configuration they hold.  It doesn't
        actually use those mechanisms because we have many config classes for
        the one `IdGenerator` class, instead of the other way around, and as a
        result a "config as factory" approach works better.
        """
        packer = self._make_dimension_packer(data_id)
        return FullIdGenerator(
            packer,
            DataCoordinate.standardize(data_id, **kwargs, graph=packer.dimensions),
            release_id=self.release_id,
            n_releases=self.n_releases,
        )

    def _make_dimension_packer(self, data_id: DataCoordinate) -> DimensionPacker:
        """Abstract hook for building a dimension packer from configuration.

        Parameters
        ----------
        data_id : `DataCoordinate`
            The data ID the `IdGenerator` will embed into all IDs.  This
            generally must be a fully-expanded data ID (i.e. have dimension
            records attached), that identifies the "instrument" or "skymap"
            dimension, though this requirement may be relaxed for certain
            dimension packer types.

        Returns
        -------
        packer : `lsst.daf.butler.DimensionPacker`
            Object that packs data IDs into integers.
        """
        raise NotImplementedError("Method is abstract.")


class DetectorExposureIdGeneratorConfig(BaseIdGeneratorConfig):
    """Configuration class for generating integer IDs from
    ``{exposure, detector}`` data IDs.

    See `IdGenerator` for usage.
    """

    packer = Instrument.make_dimension_packer_config_field()

    def _make_dimension_packer(self, data_id: DataCoordinate) -> DimensionPacker:
        # Docstring inherited.
        return self.packer.apply(data_id, is_exposure=True)


class DetectorVisitIdGeneratorConfig(BaseIdGeneratorConfig):
    """Configuration class for generating integer IDs from
    ``{visit, detector}`` data IDs.

    See `IdGenerator` for usage.
    """

    packer = Instrument.make_dimension_packer_config_field()

    def _make_dimension_packer(self, data_id: DataCoordinate) -> DimensionPacker:
        # Docstring inherited.
        return self.packer.apply(data_id, is_exposure=False)


class SkyMapIdGeneratorConfig(BaseIdGeneratorConfig):
    """Configuration class for generating integer IDs from
    ``{tract, patch, [band]}`` data IDs.

    See `IdGenerator` for usage.
    """

    packer = SkyMapDimensionPacker.make_config_field()

    def _make_dimension_packer(self, data_id: DataCoordinate) -> DimensionPacker:
        # Docstring inherited.
        return self.packer.apply(data_id)


class IdGenerator:
    """A helper class for packing some combination of a data ID, a per-data-ID
    counter, and a release ID into a single 64-bit integer.

    As an object frequently passed into code that otherwise has no knowledge of
    its own data ID, `IdGenerator` also implements ``__str__`` to provide a
    human-readable representation of the data ID for use in logs and exception
    messages, with a suitable fallback when no data ID was provided to it.

    Notes
    -----
    Instances of this class are expected to usually be created via
    configuration, which will return a derived instance.  This pattern starts
    with one of `DetectorExposureIdGeneratorConfig`,
    `DetectorVisitIdGeneratorConfig`, and `SkyMapIdGeneratorConfig` (which have
    the same interface), and looks something this:

        from lsst.meas.base import DetectorVisitIdGeneratorConfig
        from lsst.pex.config import Config
        from lsst.pipe.base import PipelineTask

        class SomeTaskConfig(PipelineTaskConfig, ...):
            id_generator = DetectorVisitIdGeneratorConfig.make_field()

        class SomeTask(PipelineTaskTask):

            ConfigClass = SomeTaskConfig

            ...

            def runQuantum(self, ..., data_id: DataCoordinate):
                id_generator = self.config.apply(data_id)
                catalog = id_generator.make_source_catalog(self.schema) ...

    There is no requirement that `IdGenerator` instances be constructed in
    `PipelineTask.runQuantum` methods and passed to the ``run`` method, but
    this is the most common approach.

    Code that wishes to instead unpack these record IDs to obtain the release
    ID, data ID and counter value should use the same config (often loaded from
    the ``Butler``) and pass a fully-expanded data ID identifying only a
    particular ``skymap`` or ``instrument`` to `unpacker_from_config`::

        config = butler.get("some_task_config")
        catalog = butler.get("some_output_catalog", given_data_id)
        unpacker = IdGenerator.unpacker_from_config(
            config.id_generator, butler.registry.expandDataId(skymap="HSC"),
        )
        release_id, embedded_data_id, counter = unpacker(catalog[0]["id"])
        assert embedded_data_id == given_data_id

    This example is a bit contrived, as the ability to reconstruct the data ID
    is really only useful when you don't have it already, such as when the
    record ID is obtained from some further-processed version of the original
    table (such as a SQL database), and in that context the right config to
    load will not be obvious unless it has been carefully documented.

    Simple instances of the base class that do not include a data ID may also
    be constructed by calling the constructor directly::

        id_generator = IdGenerator()

    These IDs may not be unpacked, but they also don't need to be, because
    they're just the per-catalog "counter" integer already.

    See Also
    --------
    :ref:`lsst.meas.base-generating-source-and-object-ids`
    """

    # TODO: remove this method on DM-38687.
    # No deprecation decorator here because the type this method accepts is
    # itself deprecated, so it's only going to be called by code paths that
    # will go away when the deprecation turns into a removal, and which already
    # warn.
    @staticmethod
    def _from_exposure_id_info(exposure_id_info: ExposureIdInfo) -> IdGenerator:
        """Construct a new ID generator from the object this class supersedes.

        This method is deprecated along with the type it accepts; it's provided
        only as a temporary helper to aid in the transition from
        `lsst.obs.base.ExposureIdInfo` to `IdGenerator`.
        """
        return _ExposureIdInfoIdGenerator(exposure_id_info)

    @property
    def catalog_id(self) -> int:
        """The integer identifier for the full catalog with this data ID, not
        just one of its rows (`int`).

        This combines the packed data ID and release ID, but not the
        counter.
        """
        return 0

    def __str__(self) -> str:
        """Return a human-readable representation of the data ID (or a note
        about its absence) for use in log and error messages.
        """
        return "[no data ID]"

    def make_table_id_factory(self) -> IdFactory:
        """Construct a new `lsst.afw.table.IdFactory` for this catalog."""
        return IdFactory.makeSimple()

    def make_source_catalog(self, schema: Schema) -> SourceCatalog:
        """Construct a empty catalog object with an ID factory.

        This is a convenience function for the common pattern of calling
        `make_table_id_factory`, constructing a `~lsst.afw.table.SourceTable`
        from that, and then constructing an (empty)
        `~lsst.afw.table.SourceCatalog` from that.
        """
        table = SourceTable.make(schema, self.make_table_id_factory())
        return SourceCatalog(table)

    def arange(self, *args, **kwargs) -> np.ndarray:
        """Generate an array of integer IDs for this catalog.

        All parameters are forwarded to `numpy.arange` to generate an array of
        per-catalog counter integers.  These are then combined with the
        `catalog_id`` to form the returned array.

        The IDs generated by `arange` will be equivalent to those generated by
        `make_table_id_factory` (and by extension, `make_source_catalog`) only
        if the counter integers start with ``1``, not ``0``, because that's
        what `~lsst.afw.table.IdFactory` does.
        """
        return np.arange(*args, **kwargs)

    @classmethod
    def unpacker_from_config(
        cls,
        config: BaseIdGeneratorConfig,
        fixed: DataCoordinate,
    ) -> Callable[[int], tuple[DataCoordinate, int]]:
        """Return a callable that unpacks the IDs generated by this class,
        from a config field.

        Parameters
        ----------
        config : `BaseIdGeneratorConfig`
            Configuration for an ID generator.
        fixed : `DataCoordinate`
            Data ID identifying the dimensions that are considered fixed by the
            `IdGenerator` that produced the IDs: usually just ``instrument`` or
            ``skymap``, depending on the configuration. For most configurations
            this will need to be a fully-expanded data ID.

        Returns
        -------
        unpacker
            Callable that takes a single `int` argument (an ID generated by an
            identically-configured `IdGenerator`) and returns a tuple of:

            - release_id: the integer that identifies a data release or
              similar (`int`);
            - data_id : the data ID used to initialize the original ID
              generator (`DataCoordinate`);
            - counter : the counter part of the original ID (`int`).

        Notes
        -----
        This method cannot be used on IDs generated without a data ID.
        """
        packer = config._make_dimension_packer(fixed)
        return cls.unpacker_from_dimension_packer(packer, config.n_releases)

    @classmethod
    def unpacker_from_dimension_packer(
        cls,
        dimension_packer: DimensionPacker,
        n_releases: int = DEFAULT_N_RELEASES,
    ) -> Callable[[int], tuple[int, DataCoordinate, int]]:
        """Return a callable that unpacks the IDs generated by this class,
        from a `lsst.daf.butler.DimensionPacker` instance.

        Parameters
        ----------
        dimension_packer : `lsst.daf.butler.DimensionPacker`
            Dimension packer used to construct the original
            `DimensionPackerIdGenerator`.
        n_releases : `int`, optional
            Number of (contiguous, starting from zero) ``release_id`` values to
            reserve space for. One (not zero) is used to reserve no space.

        Returns
        -------
        unpacker
            Callable that takes a single `int` argument (an ID generated by an
            identically-constructed `DimensionPackerIdGenerator`) and returns a
            tuple of:

            - release_id: the integer that identifies a data release or
              similar (`int`);
            - data_id : the data ID used to initialize the original ID
              generator (`DataCoordinate`);
            - counter : the counter part of the original ID (`int`).

        Notes
        -----
        This method cannot be used on IDs generated with no data ID.
        """
        bits = _IdGeneratorBits(dimension_packer, n_releases)

        def unpack(record_id: int) -> tuple[int, DataCoordinate, int]:
            rest, counter = divmod(record_id, bits.n_counters)
            rest, packed_data_id = divmod(rest, bits.n_data_ids)
            rest, release_id = divmod(rest, bits.n_data_ids)
            if rest:
                raise ValueError(
                    f"Unexpected overall factor {rest} in record_id {record_id}, "
                    f"after extracting packed_data_id={packed_data_id}, counter={counter}, and "
                    f"release_id={release_id}."
                )
            data_id = bits.packer.unpack(packed_data_id)
            return release_id, data_id, counter

        return unpack


class FullIdGenerator(IdGenerator):
    """The subclass of `IdGenerator` that actually includes packed data IDs
    and release IDs in its generated IDs.

    Parameters
    ----------
    dimension_packer : `lsst.daf.butler.DimensionPacker`
        Object that packs data IDs into integers.
    data_id : `lsst.daf.butler.DataCoordinate`
        Data ID to embed in all generated IDs and random seeds.
    release_id : `int`, optional
        Release identifier to embed in generated IDs.
    n_releases : `int`, optional
        Number of (contiguous, starting from zero) `release_id` values to
        reserve space for. One (not zero) is used to reserve no space.

    Notes
    -----
    Instances of this class should usually be constructed via configuration
    instead of by calling the constructor directly; see `IdGenerator` for
    details.
    """

    def __init__(
        self,
        dimension_packer: DimensionPacker,
        data_id: DataCoordinate,
        release_id: int = DEFAULT_RELEASE_ID,
        n_releases: int = DEFAULT_N_RELEASES,
    ):
        self._bits = _IdGeneratorBits(dimension_packer, n_releases)
        self._release_id = release_id
        self._data_id = data_id.subset(self._bits.packer.dimensions)
        self._packed_data_id = self._bits.packer.pack(self._data_id)

    @property
    def data_id(self) -> DataCoordinate:
        """The data ID that will be embedded in all generated IDs
        (`DataCoordinate`)."""
        return self._data_id

    @property
    def release_id(self) -> int:
        """The release ID that will embedded in all generated IDs (`int`)."""
        return self._release_id

    @property
    def catalog_id(self) -> int:
        # Docstring inherited.
        return self._packed_data_id + self._bits.n_data_ids * self._release_id

    def __str__(self) -> str:
        # Docstring inherited.
        return str(self.data_id)

    def make_table_id_factory(self) -> IdFactory:
        # Docstring inherited.
        return IdFactory.makeSource(self.catalog_id, self._bits.counter_bits)

    def arange(self, *args, **kwargs) -> np.ndarray:
        # Docstring inherited.
        lower = super().arange(*args, **kwargs)
        if np.any(lower >= self._bits.n_counters):
            arg_terms = [repr(arg) for arg in args] + [f"{k}={v!r}" for k, v in kwargs.items()]
            raise ValueError(
                f"Integer range from numpy.arange({arg_terms}) has "
                f"{(lower >= self._bits.n_counters).sum()} values that are not "
                f"below the upper bound of {self._bits.n_counters}."
            )
        return lower + self.catalog_id * self._bits.n_counters


@dataclasses.dataclass
class _IdGeneratorBits:
    """A private helper struct that manages the allocation of bits between the
    packed data ID, the release ID, and a per-catalog counter.
    """

    packer: DimensionPacker
    """Object that maps data IDs to integers
    (`lsst.daf.butler.DimensionPacker`).
    """

    n_releases: int = dataclasses.field(default=0)
    """Number of releases to reserve space for, starting from zero (`int`)."""

    n_data_ids: int = dataclasses.field(init=False)
    """Number of contiguous packed data IDs to reserve space for, starting
    from zero (`int`).
    """

    counter_bits: int = dataclasses.field(init=False)
    """Number of bits allocated to the per-catalog counter (`int`)."""

    n_counters: int = dataclasses.field(init=False)
    """Number of contiguous counter values to reserve space for, starting from
    zero (`int`)."""

    def __post_init__(self) -> None:
        self.n_data_ids = 1 << self.packer.maxBits
        upper_bits = (self.n_releases - 1).bit_length() + self.packer.maxBits
        self.counter_bits = IdFactory.computeReservedFromMaxBits(upper_bits)
        self.n_counters = 1 << self.counter_bits


# TODO: remove this method on DM-38687.
# No deprecation decorator here because the type this class holds is itself
# deprecated, so it's only going to be called by code paths that will go away
# when the deprecation turns into a removal, and which already warn.
class _ExposureIdInfoIdGenerator(IdGenerator):
    """A `IdGenerator` implementation to aid in the transition from
    `lsst.obs.base.ExposureIdInfo`.
    """

    def __init__(self, exposure_id_info: ExposureIdInfo):
        self._exposure_id_info = exposure_id_info

    @property
    def catalog_id(self) -> int:
        # Docstring inherited.
        return self._exposure_id_info.expId

    def __str__(self) -> str:
        return str(self.catalog_id)

    def make_table_id_factory(self) -> IdFactory:
        # Docstring inherited.
        return self._exposure_id_info.makeSourceIdFactory()

    def arange(self, *args, **kwargs) -> np.ndarray:
        # Docstring inherited.
        raise NotImplementedError(
            "This IdGenerator implementation does not support arange; "
            "please update to IdGenerator.from_config for a full-featured implementation."
        )
