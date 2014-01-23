#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import numpy

import lsst.afw.image
import lsst.afw.table
import lsst.afw.display

import lsst.meas.base

DATA_DIR = os.path.join(os.path.split(__file__)[0], "..", "tests", "data")

def main():
    calexp = lsst.afw.image.ExposureF(os.path.join(DATA_DIR, "calexp-01.fits"))
    catalog = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, "truthcat-01.fits"))
    footprints = {record.getId(): (record.getParent(), record.getFootprint())
                  for record in catalog}
    lsst.afw.display.ds9.mtv(calexp, frame=0)
    lsst.afw.display.ds9.setMaskTransparency(70)
    noiseReplacer = lsst.meas.base.NoiseReplacer(calexp, footprints, noiseSource="measure", noiseOffset=0.0,
                                                 noiseSeed=100)
    for record in catalog:
        noiseReplacer.insertSource(record.getId())
        lsst.afw.display.ds9.mtv(calexp, frame=record.getId())
        noiseReplacer.removeSource(record.getId())
    noiseReplacer.end()

if __name__ == "__main__":
    main()
