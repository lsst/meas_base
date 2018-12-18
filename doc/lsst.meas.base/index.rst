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

Task reference
==============

Command-line tasks
------------------

.. lsst-cmdlinetasks::
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

.. .. The :skip: options below exclude objects which are imported with
..    lsst.meas.base but aren't actually part of it (e.g. lsst.meas.base.Calib is
..    “actually” lsst.afw.image.calib.Calib). This causes the Sphinx build to
..    fail, at least when executed standalone.
..

.. .. automodapi:: lsst.meas.base
..    :skip: Calib
..    :skip: FatalAlgorithmError
..    :skip: MeasurementError
..    :skip: PixelValueError
..    :skip: SkyWcs

.. automodapi:: lsst.meas.base.apCorrRegistry
.. automodapi:: lsst.meas.base.applyApCorr
.. automodapi:: lsst.meas.base.baseMeasurement
.. automodapi:: lsst.meas.base.catalogCalculation
.. automodapi:: lsst.meas.base.classification
.. automodapi:: lsst.meas.base.footprintArea
.. automodapi:: lsst.meas.base.forcedMeasurement
.. automodapi:: lsst.meas.base.forcedPhotCcd
.. automodapi:: lsst.meas.base.forcedPhotImage
.. automodapi:: lsst.meas.base.measurementInvestigationLib
.. automodapi:: lsst.meas.base.noiseReplacer
.. automodapi:: lsst.meas.base.pluginRegistry

.. The skips below are necessary to prevent a build failure when building
   module-by-module, but not when building everything at once.

.. automodapi:: lsst.meas.base.plugins
   :skip: SingleFrameVariancePlugin
   :skip: ForcedVariancePlugin
   :skip: SingleFrameInputCountPlugin
   :skip: ForcedInputCountPlugin

.. automodapi:: lsst.meas.base.pluginsBase
.. automodapi:: lsst.meas.base.references
.. automodapi:: lsst.meas.base.tests
