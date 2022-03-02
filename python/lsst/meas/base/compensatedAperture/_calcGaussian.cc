#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include <math.h>
#include <vector>

namespace py = pybind11;

namespace lsst {
namespace meas {
namespace base {

py::tuple gaussian_filt_inner_product(py::array_t<double, py::array::c_style | py::array::forcecast> & array,
                                   double x_mean, double y_mean, double sig, double t) {
    //Verify the input array is conditioned in an appropriate manner.
    py::buffer_info buffer = array.request();

    if (buffer.ndim != 2) {
        throw std::runtime_error("Number of dimensions for array must be 2");
    }
    if (buffer.shape[0] != buffer.shape[1]){
        throw std::runtime_error("Array must be square");
    }
    if (!(buffer.shape[0] % 2)) {
        throw std::runtime_error("Number of pixels along array side must be odd");
    }

    double sig_sq = sig*sig;
    double t_sq = t*t;

    // fast access to array since we know all the bounds are appropriate
    auto array_unchecked = array.unchecked<2>();

    // declare most variables that will be used
    double result = 0;
    double x_offset_sq;
    double y_offset_sq;
    double y_component;
    double y_component_out;
    double tmp_x_div;
    double tmp_y_div;
    double two_sig_sq = 2*sig_sq;

    int half_domain = floor(buffer.shape[0] / 2);
    int stop = buffer.shape[0];

    double inner_norm = 1/(sig*sqrt(2*M_PI));
    double outer_norm = inner_norm/t;

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
    double sum_weights = 0;
    double weighted_sum = 0;
    double weight_square_sum = 0;

    // calculate the x profile
    for (int j = 0; j < stop; ++j) {
        double x_pos = j - half_domain;
        x_offset_sq = -1*pow(x_pos - x_mean, 2);
        tmp_x_div = x_offset_sq / two_sig_sq;
        x_container.push_back((1 / (sig * sqrt(2 * M_PI))) * exp(tmp_x_div));
        x_container_out.push_back((1 / (sig * t * sqrt(2 * M_PI))) * exp(tmp_x_div / t_sq));
    }

    // for each y row, go over the saved x vector accumulating the inner
    // product.
    for (int i = 0; i < stop; ++i) {
        double y_pos = i - half_domain;
        y_offset_sq = -1*pow(y_pos - y_mean, 2);
        tmp_y_div = y_offset_sq / two_sig_sq;
        y_component = (1 / (sig * sqrt(2 * M_PI))) * exp(tmp_y_div);
        y_component_out = (1 / (sig * t * sqrt(2 * M_PI))) * exp(tmp_y_div / t_sq);
        for (int j = 0; j < stop; ++j){
            double weight = 4*M_PI*((t_sq+1)/(t_sq-1))*
                (y_component*x_container[j] - y_component_out*x_container_out[j]);
            double weighted_value = weight*array_unchecked(i, j);
            result += weighted_value;
            if (weight < 0)  {
                // this block is used to calculate weighted variance
                weight *= -1;
                sum_weights += weight;
                weighted_sum += -1 * weighted_value;
                weight_square_sum += -1* weighted_value * array_unchecked(i, j);
            }
        }
    }
    double variance = weight_square_sum/sum_weights - (weighted_sum*weighted_sum)/(sum_weights*sum_weights);
    //return result;
    return py::make_tuple(result, variance);
}

PYBIND11_MODULE(_calcGaussian, mod) {
    mod.def("_gaussianFiltInnerProduct", &gaussian_filt_inner_product,
R"doc_string(Calculates compensated aperture for a Gaussian filter.

Parameters
----------
array : numpy array of `float`
    Array in which the Gaussian filter will be applied, must be square and
    odd shaped.
x_mean : `float`
    The x center position of the object in the reference frame of the array
    argument
y_mean : `float`
    The y center position of the object in the reference frame of the array
    argument
sig : `float`
    The std of the Gaussian kernel
t : `float`
    The dilating parameter used to calculate background.

Returns
-------
result : `float`
    Compensated aperture flux

Raises
------
RuntimeError
    Raised if the input array is not properly shaped

)doc_string");
}

} // namespace base
} // namespace meas
} // namespace lsst