#pragma once

namespace math_utils {

constexpr int min_to_m(int m) {
    return m % 2 == 0 ? 1 : -1;
}

double omega(int i, int j);
}
