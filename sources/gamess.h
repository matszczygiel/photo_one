#ifndef GAMESS_H
#define GAMESS_H

#include <complex>
#include <stdexcept>
#include <vector>

#include "constants.h"

namespace Gamess {

constexpr int crt_siz[11] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55};

int pos_change_gamess(const int &l, const int &pos);

template <class type>
std::vector<type> order(const std::vector<std::vector<type>> &shl_crt,
                        const int &l) {
    /* transform to the Gamess shell indexing */
    int p_gam, pos;
    std::vector<type> shl_crt_dum(crt_siz[l]);
    for (int p = 0; p <= l; p++)
        for (int q = 0; q <= l - p; q++) {
            pos                = (l + 1) * p - p * (p - 1) / 2 + q;
            p_gam              = pos_change_gamess(l, pos);
            shl_crt_dum[p_gam] = shl_crt[p][q];
        }

    return shl_crt_dum;
}

template <class type1, class type2>
std::vector<type1> order_set(const std::vector<std::vector<std::vector<type2>>> &shl_crt) {
    std::vector<type1> shl_crt_dum;
    int l_max = shl_crt.size() - 1;

    for (int l = 0; l <= l_max; l++) {
        std::vector<type2> mem = order(shl_crt.at(l), l);

        for (int i = 0; i < crt_siz[l]; i++)
            shl_crt_dum.push_back(mem[i]);
    }
    return shl_crt_dum;
}
};  // namespace Gamess

#endif