###############################
Investigating LSST Measurements
###############################

Before the LSST software stack can make measurements on a source identified
within an image, it must first mask out any surrounding sources so that they
do not contaminate the measurements. The software accomplishes this by
replacing every other source in an image with noise that corresponds to the
measured background distribution.

When attempting to understand the output of measurements after the analysis is
complete, it is often useful to put an exposure back into the same state it was
in when the measurement was made. What this means in practice is that the noise
used to mask out adjacent sources will need to be reconstructed in the same
manner (using the same seed) as in the original measurement. This returns the
image to a bitwise identical state, so any measurements will produce the exact
same output. The example below demonstrates how to do this using a fully
processed ci-hsc dataset.

.. code-block:: python

    from lsst.meas.base.measurementInvestigationLib import rebuildNoiseReplacer
    from lsst.daf.butler import Butler

    ciHscDataPath = "" # Set this to the path to a ci-hsc data repository.
    ciHscDataPath = "/ssd/nlust/repos_lsst/ci_hsc/DATA/rerun/ci_hsc/"

    # Create a butler object for loading in the data.
    butler = Butler(ciHscDataPath)

    # Create a data Id for a single ccd.
    dataId = {"visit": 903334, "detector": 16, "physical_filter": "HSC-R"}

    # Load in the calibrated exposure, and the associated source catalog.
    exposure = butler.get("calexp", dataId=dataId)
    srcCat = butler.get("src", dataId=dataId)

    # Reconstruct a noise replacer from the loaded data.
    noiseReplacer = rebuildNoiseReplacer(exposure, srcCat)

The  NoiseReplacer object, when created, records the pixel values for all
sources in the exposure internally and then sets the pixel values in the
exposure to random noise drawn from the background distribution. This allows
individual sources to be visualized by adding them source back in into the
exposure with a call to the NoiseReplacer's ``insertSource`` method. The source
can be  removed again with a call to the ``removeSource`` method. When all
operations involving the NoiseReplacer are completed, calling ``end`` will
reset the exposure back to its original state.

Rather than inserting and viewing individual sources, a NoiseReplacer object
is more commonly used in the context of making measurements. The
following example shows how to rerun measurements, using a reconstructed
NoiseReplacer. The example selects a subset of objects for which measurements
are to be rerun, and places the outputs into a new catalog, retaining the same
id for each object.

.. code-block:: python

    # Continued from the above example
    from lsst.afw.table import SourceTable
    from lsst.meas.base.measurementInvestigationLib import makeRerunCatalog
    from lsst.meas.base import (SingleFrameMeasurementConfig,
                                SingleFrameMeasurementTask)

    # Make a list of ids of objects to remeasure
    idsToRerun = [775958066192449538, 775958066192449539,
                  775958066192449540, 775958066192449541]

    # Fields to copy from old catalog, these are generally fields added outside
    # the measurement framework, that may be desirable to maintain
    fields = ["deblend_nChild"]

    # Create a new schema object, and use it to initialize a measurement task
    schema = SourceTable.makeMinimalSchema()

    # Configure any plugins at this stage.
    measConfig = SingleFrameMeasurementConfig()

    measTask = SingleFrameMeasurementTask(schema, config=measConfig)

    # Create a Measurement catalog containing only the ids to remeasure
    newSrcCatalog = makeRerunCatalog(schema, srcCat, idsToRerun, fields=fields)

    # Re-run measure on the sources selected above, using the reconstructed
    # noise replacer.
    measTask.runPlugins(rebuildNoiseReplacer(exposure, srcCat), newSrcCatalog, exposure)


At the time of writing, measTask.runPlugins will not run on any child objects
if their parents are not also in the catalog. The makeRerunCatalog function
has two options to resolve this issue. The ``resetParents`` flag (``True`` by
default) will reset the parent keys of any children whose parents are not
included in the idList, effectively turning them into parents. The
``addParents`` flag (False by default) will add these parents to the idList.

.. code-block:: python

    import numpy as np
    # Get only the child ids
    idsToRerun = srcCat["id"][srcCat[srcCat.getParentKey()] != 0]

    # Previously, it would skip all children - only the noiseReplacer takes time
    newSrcCatalog = makeRerunCatalog(
        schema, srcCat, idsToRerun, fields=fields, resetParents=False
    )
    measTask.runPlugins(rebuildNoiseReplacer(exposure, srcCat), newSrcCatalog, exposure)
    # None of these objects have centroids
    print(len(newSrcCatalog), np.sum(np.isnan(newSrcCatalog["slot_Centroid_x"])))

    # resetParents=True (default) resets parents and takes a few seconds longer
    newSrcCatalog = makeRerunCatalog(
        schema, srcCat, idsToRerun, fields=fields, resetParents=True
    )
    measTask.runPlugins(rebuildNoiseReplacer(exposure, srcCat), newSrcCatalog, exposure)
    # Now none of the objects have nan centroids
    print(len(newSrcCatalog), np.sum(np.isnan(newSrcCatalog["slot_Centroid_x"])))

    # Setting addParents=True adds all parents and takes a little longer still
    newSrcCatalog = makeRerunCatalog(
        schema, srcCat, idsToRerun, fields=fields, addParents=True
    )
    measTask.runPlugins(rebuildNoiseReplacer(exposure, srcCat), newSrcCatalog, exposure)
    # None of the objects have nan centroids and the catalog is larger than above
    print(len(newSrcCatalog), np.sum(np.isnan(newSrcCatalog["slot_Centroid_x"])))
