// -*- lsst-c++ -*-

#include "boost/make_shared.hpp"

#include "lsst/meas/base/PixelFlags.h"
#include "lsst/afw/detection/FootprintFunctor.h"

namespace lsst { namespace meas { namespace base { namespace algorithms {

template <typename MaskedImageT>
class FootprintBits : public afw::detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintBits(MaskedImageT const& mimage) :
        afw::detection::FootprintFunctor<MaskedImageT>(mimage), _bits(0)
    {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _bits = 0x0;
    }
    virtual void reset(afw::detection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        _bits |= loc.mask(0, 0);
    }

    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    typename MaskedImageT::Mask::Pixel _bits;
};


}}}} // namespace lsst::meas::base::algorithms
