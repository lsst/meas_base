from lsst.meas.base.sfm import SingleFrameMeasurementTask
root.tableVersion=0
root.measurement.algorithms = ["shape.sdss", "flux.gaussian", "flux.psf", "flux.sinc", "flux.naive", "skycoord", "classification.extendedness","flags.pixel"]
root.measurement.slots.centroid = "centroid.sdss"
root.measurement.slots.shape = "shape.sdss"
root.measurement.slots.apFlux = "flux.sinc"
root.measurement.slots.modelFlux = "flux.gaussian"
root.measurement.slots.psfFlux = "flux.psf"
root.measurement.slots.instFlux = "flux.naive"
