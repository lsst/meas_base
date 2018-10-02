##################################
Overview of the measurement system
##################################

The meas_base package provides the basis of LSST's system for performing measurements on astronomical images.
This document provides a brief overview of how the pieces of the meas_base framework fit together.
It is intended to help new users orient themselves rather than to provide a comprehensive user guide.

The structure of meas_base measurement
======================================

Measurements take place within the meas_base framework.

meas_base distinguishes between a measurement “algorithm” (i.e. a bit of code which takes some pixels and returns some useful result based upon them) and a measurement “driver”, which takes an image and runs a bunch of algorithms on it.

Typical measurement algorithms include (say) PSF fluxes, the SDSS shape algorithm, etc.

The measurement driver depends on the type of measurement being performed.
meas_base provides two: single frame measurement, which is measurement on a single visit, and forced measurement, which takes an image and an external catalog which it uses to fix source centroids and shapes.

Measurement drivers are implemented as tasks, while algorithms are meas_base plugins.
In most (but not all) cases, algorithms and drivers can be compiled fairly arbitrarily; that is, the PSF flux algorithm can be used with single frame measurement or forced measurement, etc.

meas_base provides those two drivers and a bunch of simple algorithms.
Other packages can provide more of either, following the same conventions.
Most of our meas_extensions packages provide measurement algorithms for use in the meas_base framework (e.g. simpleShape, shapeHSM, etc).

A notable extension package is ip_diffim.
This provides a driver suitable for meas_base-style measurement on difference images and a bunch of algorithms (e.g. dipole measurement) which are of relevance in that context.

Naming and selecting measurement algorithms
===========================================

The new users can generally skip over the details of measurement drivers: we have the ones we need for now, and they mostly “just work”.
The interesting stuff all happens in algorithms.

Understand first that algorithms may be referred to by strings which describe both their home package and their name.
Thus, ``base_PsfFlux`` refers to the ``PsfFlux`` algorithm in :py:mod:`lsst.meas.base`; ``ip_diffim_DipoleFit`` refers to the ``DipoleFit`` algorithm in :py:mod:`lsst.ip.diffim`.

When configuring a measurement driver (i.e. a task like :py:class:`lsst.meas.base.sfm.SingleFrameMeasurementTask`) we specify a list of plugins to be run when executing the driver.
See, for example, `here`_.
These can be over-ridden by the user in the task configuration.

When the driver runs, it will return an ``afw::table``.
Each row in the table contains the results of running all of the plugins which have been enabled on a particular source.
Each algorithm plugin can add columns to the tables to contain its results.
These columns build upon the naming scheme described above.
Thus, if the ``base_PsfFlux`` algorithm is run, you will see a bunch of fields in the table with names starting ``base_PsfFlux_`` — ``base_PsfFlux_instFlux``, ``base_PsfFlux_instFluxErr``, and so on.

(Incidentally, note that this means that the schema of the output table depends on the particular plugins that were enabled at task configuration time. This makes it awkward to provide a single document that describes pipeline outputs — they vary with user configuration.)

Algorithm fields and flags
==========================

Following from the above, when a new measurement driver (task) is created, it will:

#. Create a new schema for its output;
#. Hand that to each of the measurement algorithms (plugins) it's planning to
   run, and ask them to add fields which will be used to store their outputs.

In the limit of poor documentation, we can look at the code for the algorithm to understand what fields it's adding and what they mean.
(Alternatively, we can inspect the schema of the output table.)

Note that algorithms can be implemented in C++ or Python, so we need to be able to look at both.

`Here is`_ an example in Python: the ``base_FPPosition`` plugin.
When this plugin is created, it is given a ``schema`` argument, which will ultimately be the schema of the measurement driver outputs.
It adds three things to the schema: fields to contain its outputs, and two flag fields::

   schema.addField(name + "_flag", type="Flag", doc="Set to True for any fatal failure")
   schema.addField(name + "_missingDetector_flag", type="Flag", doc="Set to True if detector object is missing")

When you look at the output of any measurement driver that has executed ``base_FPPosition``, you should see those fields in the schema.
And you can see that the flags are defined and (briefly) documented: if the algorithm fails for any reason, ``base_FPPosition_flag`` will be set to ``True``, and if it fails specifically because the detector is missing, ``base_FPPosition_missingDetector_flag`` will be set to ``True``.

Several other plugins are defined in Python and follow the same pattern.  See, for example, ``base_Jacobian`` (`look here <https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/plugins.py#L156>`_), ``base_PeakCentroid`` (`and here <https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/plugins.py#L319>`_), etc.

Of course, plugins which are complex or have to run fast are written in C++.
These follow the same general principles, except that the code may be a little less familiar looking.
Check out the flags defined for `aperture`_ or `PSF`_ fluxes, for example.

Flags based on masks
====================

Often we don't care about failures with a particular measurement algorithm, which will be flagged as above, but about whether a source happened to fall on some area of the image which we know is suspect (interpolated, at the edge of the CCD, saturated, whatever).

This information is captured at the pixel level by the image mask planes.
A special measurement plugin, ``base_PixelFlags``, can be run to transfer information from the image mask to the source flags.
For example, when the ``base_PixelFlags`` algorithm is run, it will set a ``base_PixelFlags_flag_interpolated`` field to ``True`` if any of the pixels in the footprint of the source being measured were marked as ``INTRP`` on the image.

.. _here: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/sfm.py#L121
.. _Here is: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/python/lsst/meas/base/plugins.py#L130
.. _aperture: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/src/ApertureFlux.cc#L44
.. _PSF: https://github.com/lsst/meas_base/blob/35d32cdfa0559496d21b7de0310bd9161e120578/src/PsfFlux.cc#L44
