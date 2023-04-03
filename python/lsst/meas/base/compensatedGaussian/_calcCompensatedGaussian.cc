#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "lsst/cpputils/python.h"
#include <math.h>
#include <vector>
#include <stdio.h>

namespace py = pybind11;

namespace lsst {
namespace meas {
namespace base {

py::tuple compensated_gaussian_filt_inner_product(py::array_t<double, py::array::c_style | py::array::forcecast> & array,
                                                  py::array_t<double, py::array::c_style | py::array::forcecast> & variance_array,
                                                  double x_mean, double y_mean, double sig, double t) {
    //Verify the input array is conditioned in an appropriate manner.
    py::buffer_info buffer = array.request();
    py::buffer_info variance_buffer = variance_array.request();

    if (buffer.ndim != 2) {
        throw std::runtime_error("Number of dimensions for array must be 2");
    }
    if (buffer.shape[0] != buffer.shape[1]){
        throw std::runtime_error("Array must be square");
    }
    if (!(buffer.shape[0] % 2)) {
        throw std::runtime_error("Number of pixels along array side must be odd");
    }

    if (variance_buffer.ndim != buffer.ndim) {
        throw std::runtime_error("Variance array must have same dimensions as image array");
    }
    if (variance_buffer.shape[0] != buffer.shape[0] || variance_buffer.shape[1] != buffer.shape[1]) {
        throw std::runtime_error("Variance array must have same dimensions as image array");
    }

    double sig_sq = sig*sig;
    double t_sq_inv = 1./(t*t);
    double t_inv = 1./t;

    // fast access to array since we know all the bounds are appropriate
    auto array_unchecked = array.unchecked<2>();
    auto variance_array_unchecked = variance_array.unchecked<2>();

    // declare most variables that will be used
    double flux = 0;
    double x_offset_sq;
    double y_offset_sq;
    double y_component;
    double y_component_out;
    double tmp_x_div;
    double tmp_y_div;
    double two_sig_sq = 2*sig_sq;

    int half_domain = floor(buffer.shape[0] / 2);
    int stop = buffer.shape[0];

    // adjust the x and y center to be centered about the middle of the array
    x_mean -= (double) half_domain;
    y_mean -= (double) half_domain;

    /*
    create containers to store the 1 dimensional x calculations, these will be
    accessed for each loop in y, no need to do the same calculations each time
    */
    std::vector<float> x_container;
    std::vector<float> x_container_out;
    x_container.reserve(stop);
    x_container_out.reserve(stop);

    // for weighted variance calculation
    double variance = 0;
    double var_norm = 0;

    // calculate the x profile
    for (int j = 0; j < stop; ++j) {
        double x_pos = j - half_domain;
        x_offset_sq = -1*pow(x_pos - x_mean, 2);
        tmp_x_div = x_offset_sq / two_sig_sq;
        x_container.push_back(exp(tmp_x_div));
        x_container_out.push_back(t_inv * exp(tmp_x_div * t_sq_inv));
    }

    // for each y row, go over the saved x vector accumulating the inner
    // product.
    for (int i = 0; i < stop; ++i) {
        double y_pos = i - half_domain;
        y_offset_sq = -1*pow(y_pos - y_mean, 2);
        tmp_y_div = y_offset_sq / two_sig_sq;
        y_component = exp(tmp_y_div);
        y_component_out = t_inv * exp(tmp_y_div * t_sq_inv);
        for (int j = 0; j < stop; ++j){
            double weight = (y_component*x_container[j] - y_component_out*x_container_out[j]);
            double weighted_value = weight*array_unchecked(i, j);
            flux += weighted_value;

            variance += weight*weight*variance_array_unchecked(i, j);
            var_norm += weight*weight;
        }
    }

    // Normalization of the normalized Gaussian filter is 4 * pi * sig**2 * (t**2 + 1)/(t**2 - 1)
    // We have deliberately not applied the normalization factor 1. / (2 * pi * sig**2) in the
    // inner product, leaving the normalization term as 2 * (t**2 + 1)/(t**2 - 1) for the
    // flux, and XXX for the variance.
    flux *= 2 * (t*t + 1) / (t*t - 1);

    // The variance has been normalized and so requires the full normalization term.
    variance /= var_norm;
    variance *= (4 * M_PI * sig_sq * (t*t + 1) ) / (t*t);

    return py::make_tuple(flux, variance);
}

void wrapCalcCompensatedGaussian(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrap([](auto &mod) {
        mod.def("_compensatedGaussianFiltInnerProduct", &compensated_gaussian_filt_inner_product,
R"doc_string(Calculates flux and variance for a Compensated Gaussian filter.

Parameters
----------
array : `np.ndarray` (N, N)
    Array in which the Gaussian filter will be applied, must be square and
    odd shaped.
variance_array : `np.ndarray` (N, N)
    Variance array in which the Gaussian filter will be applied, must be
    same shape as array.
x_mean : `float`
    The x center position of the object in the reference frame of the array
    argument.
y_mean : `float`
    The y center position of the object in the reference frame of the array
    argument.
sig : `float`
    The std of the Gaussian kernel
t : `float`
    The dilating parameter used to calculate background.

Returns
-------
result : `tuple` of `float`
    Compensated Gaussian flux and variance.

Raises
------
RuntimeError
    Raised if the input array is not properly shaped

)doc_string");
    });
}

} // namespace base
} // namespace meas
} // namespace lsst
