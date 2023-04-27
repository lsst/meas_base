.. py:currentmodule:: lsst.meas.base

.. _lsst.meas.base-generating-source-and-object-ids:

################################
Generating source and object IDs
################################

The tables used to hold source and object measurements are expected to have row IDs that are unique across an entire processing run.
This is achieved by first packing the :ref:`Butler data ID <lsst.daf.butler-dimensions_data_ids>` that identifies an in-memory table object into an integer, and then packing that with a per-table counter.
The packing algorithms are reversible (they're not just a hash), and we always aim to fit the full row ID into a single 64-bit signed integer.
In special cases, like LSST data releases, an identifier for the release may be packed in as well, making the row ID globally unique (for its type) across all releases.

The `IdGenerator` class is the main entry point for generating these IDs.
Usage starts with with defining an `lsst.pex.config` field for the type of data ID that identifies a table object:

- `SkyMapIdGeneratorConfig` for ``{tract, patch}`` or ``{tract, patch, band}}`` data IDs (most object tables);
- `DetectorVisitIdGeneratorConfig` for ``{visit, detector}`` data IDs (most source tables);
- `DetectorExposureIdGeneratorConfig` for ``{exposure, detector}`` data IDs (rare, at least in the default LSST data model).

These have a `~BaseIdGeneratorConfig.make_field` method that can be used to define the config field with minimal boilerplate.

After configuration, the `~BaseIdGeneratorConfig.apply` method can then be called to make an `IdGenerator` instance.
`IdGenerator` instances can be used to make `lsst.afw.table.IdFactory` instances via `IdGenerator.make_table_id_factory` or `IdGenerator.make_source_catalog` or equivalent numpy arrays of IDs via the `IdGenerator.arange` method.
The `~IdGenerator.catalog_id` property provides access to just the packed data ID, which can be useful as a random number generator seed or an ID for the rows of summary tables whose rows correspond to images, not sources or objects.
See the `IdGenerator` class documentation for details and some examples.
