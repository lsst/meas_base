#include "lsst/afw/table.h"
#include "lsst/afw/image.h"
#include "lsst/meas/base/SincFlux.h"

using namespace lsst::afw::table;
using namespace lsst::afw::image;
using namespace lsst::meas::base;

int main() {
    auto truthCat = SourceCatalog::readFits("truth.fits");
    auto exposure = Exposure<float>("exposure.fits");
    SchemaMapper mapper(truthCat.getSchema());
    mapper.addMinimalSchema(truthCat.getSchema(), true);
    SincFluxControl ctrl;
    ctrl.radius1 = 0.0;
    ctrl.radius2 = 16.0;
    Key<double> xKey = mapper.editOutputSchema().addField<double>("peak_x", "peak column position", "pixels");
    Key<double> yKey = mapper.editOutputSchema().addField<double>("peak_y", "peak row position", "pixels");
    PTR(SingleFrameAlgorithm) algorithm(
        new SincFluxAlgorithm(ctrl, "base_SincFlux", mapper.editOutputSchema())
    );
    SourceCatalog measCat(mapper.getOutputSchema());
    measCat.getTable()->defineCentroid("peak");
    measCat.insert(mapper, measCat.begin(), truthCat.begin(), truthCat.end());
    for (auto & measRecord : measCat) {
        measRecord.set(xKey, measRecord.getFootprint()->getPeaks().front()->getFx());
        measRecord.set(yKey, measRecord.getFootprint()->getPeaks().front()->getFy());
        algorithm->measure(measRecord, exposure);
    }
    return 0;
}
