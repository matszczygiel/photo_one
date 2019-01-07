#include "harmonics.h"

#include <cmath>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "math.h"

using boost::math::factorial;
using boost::math::binomial_coefficient;

using namespace math_utils;


double Harmonics::NoNormCalcClmR(const int &l, const int &m,
                                 const int &lx, const int &ly, const int &lz) {
    if (lx + ly + lz != l)
        return 0.0;

    if (m > 0) {
        int am   = std::abs(m);
        int bs   = l - lz - am;
        int twol = 2 * l;
        if (bs % 2 == 1 || bs < 0)
            return 0.0;
        bs  = bs / 2;
        double fac = M_SQRT2 * omega(l, m) * min_to_m(am) /
                     (std::pow(2.0, l) * factorial<double>(l));

        int lima  = l - am;
        lima  = (lima % 2 == 0) ? lima / 2 : (lima - 1) / 2;
        int lim1b = lx - am;
        lim1b = (lim1b % 2 == 0) ? lim1b / 2 : (lim1b + 1) / 2;
        lim1b = std::max(lim1b, 0);
        int lim2b = (lx % 2 == 0) ? lx / 2 : (lx - 1) / 2;
        lim2b = std::min(lim2b, bs);

        double val = 0.0;
        for (int a = std::max(bs, 0); a <= lima; a++) {
            int twoa = 2 * a;
            double faci = binomial_coefficient<double>(l, a) * min_to_m(a) *
                          factorial<double>(twol - twoa) * binomial_coefficient<double>(a, bs) /
                          factorial<double>(l - am - twoa);
            for (int b = lim1b; b <= lim2b; b++) {
                int twob = 2 * b + am - lx;
                if (twob % 2 == 1)
                    continue;
                val += binomial_coefficient<double>(bs, b) *
                       binomial_coefficient<double>(am, twob) *
                       min_to_m(twob / 2) * faci;
            }
        }
        return min_to_m(am) * val * fac / std::sqrt(2.0 * M_PI);
    }
    if (m == 0) {
        int lxx = lx;
        int lxy = lx + ly;
        if (lxx % 2 == 1)
            return 0.0;
        if (lxy % 2 == 1)
            return 0.0;
        int twol = 2 * l;
        lxx /= 2;
        lxy /= 2;
        double fac = std::sqrt((2 * l + 1) / 2.) * binomial_coefficient<double>(lxy, lxx) /
                     (std::pow(2.0, l) * factorial<double>(l));
        int lim1a = lxy;
        int lim2a = (l % 2 == 0) ? l / 2 : (l - 1) / 2;

        double val = 0.0;
        for (int a = lim1a; a <= lim2a; a++) {
            int twoa = 2 * a;
            val += binomial_coefficient<double>(l, a) * min_to_m(a) *
                   factorial<double>(twol - twoa) * binomial_coefficient<double>(a, lxy) /
                   factorial<double>(l - twoa);
        }
        return val * fac / std::sqrt(2.0 * M_PI);
    }
    if (m < 0) {
        int am   = std::abs(m);
        int bs   = l - lz - am;
        int twol = 2 * l;
        if (bs % 2 == 1 || bs < 0)
            return 0.0;
        bs  = bs / 2;
        double fac = -M_SQRT2 * omega(l, am) * min_to_m(am) /
                     (std::pow(2.0, l) * factorial<double>(l));

        int lima  = l - am;
        lima  = (lima % 2 == 0) ? lima / 2 : (lima - 1) / 2;
        int lim1b = lx - am;
        lim1b = (lim1b % 2 == 0) ? lim1b / 2 : (lim1b + 1) / 2;
        lim1b = std::max(lim1b, 0);
        int lim2b = (lx % 2 == 0) ? lx / 2 : (lx - 1) / 2;
        lim2b = std::min(lim2b, bs);

        double val = 0.0;
        for (int a = std::max(bs, 0); a <= lima; a++) {
            int twoa = 2 * a;
            double faci = binomial_coefficient<double>(l, a) * min_to_m(a) *
                          factorial<double>(twol - twoa) * binomial_coefficient<double>(a, bs) /
                          factorial<double>(l - am - twoa);
            for (int b = lim1b; b <= lim2b; b++) {
                int twob = 2 * b + am - lx + 1;
                if (twob % 2 == 1)
                    continue;
                val += binomial_coefficient<double>(bs, b) *
                       binomial_coefficient<double>(am, twob - 1) *
                       min_to_m(twob / 2) * faci;
            }
        }
        return min_to_m(am) * val * fac / std::sqrt(2.0 * M_PI);
    }
    return 0;
}

double Harmonics::NoNormYrl(const int &l, const int &m, const Eigen::Vector3d &r) {
    double res = 0;

    for (int lx = 0; lx <= l; ++lx)
        for (int ly = 0; ly <= l - lx; ++ly) {
            int lz = l - lx - ly;
            res += NoNormCalcClmR(l, m, lx, ly, lz) *
                   std::pow(r(0), lx) * std::pow(r(1), ly) * std::pow(r(2), lz);
        }
    return res;
}

double Harmonics::Y(const int &l, const int &m, const Eigen::Vector3d &r) {
    return NoNormYrl(l, m, r) / std::pow(r.norm(), l);
}
