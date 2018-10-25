#####################################
Overview of the measurement framework
#####################################

.. The meas_base package provides the basis of LSST's system for performing measurements on astronomical images.
   This page provides an overview of how the pieces of the meas_base framework fit together.
   It is intended to help new users orient themselves rather than to provide a comprehensive user guide.

The structure of the measurement framework
==========================================

`lsst.meas.base` provides the basis of LSST's system for performing measurements on astronomical images.

Measurements of astronomical objects are made by applying a particular *algorithm* — for example, an algorithm which measures PSF flux, or characterizes shape as an elliptical Gaussian — to a collection of pixels representing the object
In the LSST measurement framework, “plugins” are used to implement particular algorithms.
The `lsst.meas.base` package provides a number of common plugins, and more can be provided in other packages.
In general, complex measurement algorithms are implemented as plugins in `lsst.meas.extensions` packages, such as `lsst.meas.extensions.simpleShape`, `lsst.meas.extensions.shapeHSM`, and so on.

Measurement may take place in a number of different contexts.
For example, a measurement may take place on the image on which a collection of sources have just been identified: this is “single frame measurement.”
Alternatively, measurement may be carried out on sources detected on an image *other than the one being measured*: this is “forced measurement.”
`lsst.meas.base` provides tasks which carry out measurement in both of the above contexts, and, once again, more can be provided in other packages which follow the `lsst.meas.base` conventions..

A particularly notable extension packages is `lsst.ip.diffim`.
This provides as a task suitable for carrying out measurement on difference images, as well as a number of plugins, such as dipole measurement, which are of particular relevance in that context.

In many cases, plugins and measurement tasks can be combined fairly arbitrarily; that is, the PSF flux plugin can be used with single frame measurement or forced measurement, etc.

Using measurement plugins
=========================

After selecting the appropriate measurement task, it is necessary to select which plugins it will run.
This is done using the regular task configuration system.

Understand first that plugins may be referred to by strings which describe both their home package and their name.
Thus, ``base_PsfFlux`` refers to the ``PsfFlux`` plugin in `lsst.meas.base`; ``ip_diffim_DipoleFit`` refers to the ``DipoleFit`` plugin in `lsst.ip.diffim`.

When configuring a measurement task (for example, `~lsst.meas.base.sfm.SingleFrameMeasurementTask`), use these names to specify a list of plugins to execute. in the ``plugins`` configuration field.
Refer to `the source`_ for an example.
These can be overridden by the user in the task configuration.

When the measurement task runs, it will return an instance of `lsst.afw.table.SourceTable`.
Each row in the table contains the results of running all of the plugins which have been enabled on a particular source.
Each plugin can add columns to the table to contain its results.
These columns build upon the naming scheme described above.
Thus, if the ``base_PsfFlux`` plugin is run, you will see fields in the table with names starting ``base_PsfFlux_`` — ``base_PsfFlux_instFlux``, ``base_PsfFlux_instFluxErr``, and so on.

.. note::

   Incidentally, this means that the schema of the output table depends on the particular plugins that were enabled at task configuration time.
   This makes it awkward to provide a single document that describes pipeline outputs — they vary with user configuration.

Output fields and flags
=======================

When a new measurement task is created, it will:

#. Create a new schema for its output;
#. Hand that to each of the measurement plugins it is configured to run, and ask them to add fields which will be used to store their outputs.

At present, the outputs of plugins aren't fully documented on the web.
Until they are, you can examine the code for the plugin to understand what fields it adds and what they mean.

Inspecting the code of Python plugins
--------------------------------------

`Here is`_ an example in Python: the ``base_FPPosition`` plugin.
When this plugin is created, it is given a ``schema`` argument, which will ultimately be the schema of the measurement task outputs.
It adds three things to the schema: fields to contain its outputs, and two flag fields::

   schema.addField(name + "_flag", type="Flag", doc="Set to True for any fatal failure")
   schema.addField(name + "_missingDetector_flag", type="Flag", doc="Set to True if detector object is missing")

When you look at the output of any measurement task that has executed ``base_FPPosition``, you will see those fields in the schema.
And you can see that the flags are defined and (briefly) documented: if the plugin fails for any reason, ``base_FPPosition_flag`` will be set to ``True``, and if it fails specifically because the detector is missing, ``base_FPPosition_missingDetector_flag`` will be set to ``True``.

Several other plugins are defined in Python and follow the same pattern.
See, for example, ``base_Jacobian`` (`look here <https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/plugins.py#L156>`_), ``base_PeakCentroid`` (`and here <https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/plugins.py#L319>`_), etc.

Inspecting the code of C++ plugins
----------------------------------

Of course, plugins which are complex or have to run fast are written in C++.
These follow the same general principles, except that the code may be a little less familiar looking.
Check out the flags defined for `aperture`_ or `PSF`_ fluxes, for example.

Flags based on masks
====================

Often you don't care about failures with a particular measurement plugin, which will be flagged as above, but about whether a source happened to fall on some area of the image which you know is suspect (interpolated, at the edge of the CCD, saturated, …)

This information is captured at the pixel level by the image mask planes.
A special measurement plugin, ``base_PixelFlags``, can be run to transfer information from the image mask to the source flags.
For example, when the ``base_PixelFlags`` plugin is run, it will set a ``base_PixelFlags_flag_interpolated`` field to ``True`` if any of the pixels in the footprint of the source being measured were marked as ``INTRP`` on the image.

.. _the source: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/sfm.py#L121
.. _Here is: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/plugins.py#L130
.. _aperture: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/src/ApertureFlux.cc#L44
.. _PSF: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/src/PsfFlux.cc#L44
