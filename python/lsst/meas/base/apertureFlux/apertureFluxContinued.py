from __future__ import absolute_import, division, print_function

from lsst.utils import continueClass

from .apertureFlux import ApertureFluxAlgorithm, ApertureFluxControl

@continueClass
class ApertureFluxAlgorithm:
    Control = ApertureFluxControl

