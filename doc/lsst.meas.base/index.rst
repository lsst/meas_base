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
   tasks_and_algorithms
   generating_source_and_object_ids
   reconstructing_measurements

Contributing
============

``lsst.meas.base`` is developed at https://github.com/lsst/meas_base.
You can find Jira issues for this module under the `meas_base <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20meas_base>`_ component.

Task reference
==============

Pipeline tasks
--------------

.. lsst-pipelinetasks::
   :root: lsst.meas.base

Tasks
-----

.. lsst-tasks::
   :root: lsst.meas.base
   :toctree: tasks

Configs
-------

.. lsst-configs::
   :root: lsst.meas.base
   :toctree: configs

Python API reference
====================

.. There are two reasons for the :skip: parameters below.
   - Objects which are imported with lsst.meas.base but aren't actually part
     of it --- e.g. lsst.meas.base.PhotoCalib is actually
     lsst.afw.image.photoCalib.PhotoCalib --- cause the Sphinx build to fail
     when executed standalone.
   - SingleFrameFromGenericPlugin and ForcedFromGenericPlugin cause the
     construction of an inheritance diagram to fail, per DM-15461.
     Unfortunately, including them here causes a warning, but that seems to be
     unavoidable.

.. automodapi:: lsst.meas.base
   :no-main-docstr:
   :skip: PhotoCalib
   :skip: FatalAlgorithmError
   :skip: MeasurementError
   :skip: PixelValueError
   :skip: SkyWcs
   :skip: SingleFrameFromGenericPlugin
   :skip: ForcedFromGenericPlugin
