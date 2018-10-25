.. py:currentmodule:: lsst.meas.base

.. _lsst.meas.base:

##############
lsst.meas.base
##############

``lsst.meas.base`` provides core astronomical measurement algorithms and base classes.

.. _lsst.meas.base-using:

Using lsst.meas.base
====================

.. toctree::
   :maxdepth: 1

   intro
   reconstructing-measurements

Contributing
============

``lsst.meas.base`` is developed at https://github.com/lsst/meas_base.
You can find Jira issues for this module under the `meas_base <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20meas_base>`_ component.

Python API reference
====================

.. The :skip: options below exclude objects which are imported with
   lsst.meas.base but aren't actually part of it (e.g. lsst.meas.base.Calib is
   “actually” lsst.afw.image.calib.Calib). This causes the Sphinx build to
   fail, at least when executed standalone.

.. automodapi:: lsst.meas.base
   :skip: Calib
   :skip: FatalAlgorithmError
   :skip: MeasurementError
   :skip: PixelValueError
   :skip: SkyWcs
