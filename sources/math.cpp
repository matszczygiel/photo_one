#include "math.h"

#include <boost/math/special_functions/factorials.hpp>

using boost::math::factorial;

double math_utils::omega(int i, int j) {
    return std::sqrt((2.0 * i + 1) * factorial<double>(i - j) /
                     (2.0 * factorial<double>(i + j)));
}
