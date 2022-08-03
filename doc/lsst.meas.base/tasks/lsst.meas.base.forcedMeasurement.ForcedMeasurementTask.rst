.. currentmodule:: lsst.meas.base
.. lsst-task-topic:: lsst.meas.base.forcedMeasurement.ForcedMeasurementTask

#####################
ForcedMeasurementTask
#####################

``ForcedMeasurementTask`` runs a set of (user-selectable) plugins, implementing particular measurement algorithms, to measure the properties of sources on a single image.
A reference catalog, containing measurements made elsewere, is used to constrain some aspects of the algorithms.

Processing Summary
==================

The set of plugins to run is set in the task configuration.
Each plugin defines the values that it measures (which correspond to columns in the output table), and then conducts measurement on each detected source.
See `ForcedPlugin` for details.
This task intializes the set of plugins (thereby defining the catalog schema) from its configuration, then invokes each plugin on each source.

Most of the time, `ForcedMeasurementTask` will be used via one of the subclasses of :lsst-task:`lsst.meas.base.forcedPhotCcd.ForcedPhotCcdTask` or :lsst-task:`lsst.meas.base.forcedPhotCoadd.ForcedPhotCoaddTask`.
These combine this measurement subtask with a ``references`` subtask (see ``lsst.meas.base.references.BaseReferencesTask`` and ``~lsst.meas.base.references.CoaddSrcReferencesTask``) to perform forced measurement using measurements performed on another image as the references.
There is generally little reason to use `ForcedMeasurementTask` outside of one of these drivers, unless it is necessary to avoid using the Butler for I/O.

Forced measurement means that the plugins are provided with a reference source containing centroid and/or shape measurements that they may use however they see fit.
Some plugins can use these to set the location and size of apertures, but others may choose to ignore this information, essentially performing an unforced measurement starting at the position of the reference source (which may nevertheless be useful for certain
investigations).
Knowing how the plugin uses the reference information is essential to interpreting its resulting measurements.
Typically, centroid and shape measurement plugins (e.g., ``SdssCentroid`` and ``SdssShape``) are performing unforced measurements.

Python API summary
==================

.. lsst-task-api-summary:: lsst.meas.base.forcedMeasurement.ForcedMeasurementTask

Configuration fields
====================

.. lsst-task-config-fields:: lsst.meas.base.forcedMeasurement.ForcedMeasurementTask
