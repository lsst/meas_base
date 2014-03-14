#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/base/Results.h"

#ifndef LSST_MEAS_BASE_ApFluxComponent_h_INCLUDED
#define LSST_MEAS_BASE_ApFluxComponent_h_INCLUDED

namespace lsst { namespace meas { namespace base {

struct ApFluxComponent {
    Flux flux; ///< Measured flux in DN.
    ErrElement fluxSigma; ///< 1-Sigma error (sqrt of variance) on flux in DN.

    ApFluxComponent();
};

}}};

#endif // LSST_MEAS_BASE_ApFluxComponent_h_INCLUDED
