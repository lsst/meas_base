from lsst.meas.base.sfm import SingleFrameMeasurementTask
root.measurement.retarget(SingleFrameMeasurementTask)
root.measurement.plugins = ["centroid.peak", "base_SdssCentroid", "base_SdssShape", "base_GaussianFlux", "base_PsfFlux", "base_SincFlux", "base_NaiveFlux", "skycoord", "classification", "base_PixelFlags"]
root.tableVersion=1
root.measurement.slots.centroid = "base_SdssCentroid"
root.measurement.slots.shape = "base_SdssShape"
root.measurement.slots.apFlux = "base_SincFlux"
root.measurement.slots.modelFlux = "base_GaussianFlux"
root.measurement.slots.psfFlux = "base_PsfFlux"
root.measurement.slots.instFlux = "base_NaiveFlux"
