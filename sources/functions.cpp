#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/math/special_functions/factorials.hpp>

#include "constants.h"
#include "functions.h"
#include "harmonics.h"

using namespace std;

double dsigma(const double& photon, const Eigen::Vector3d& polarization, const Eigen::Vector3cd& dipole) {
    double pre_fct  = 4. * M_PI * M_PI * photon / speed_of_light;
    double post_fct = std::norm(polarization.dot(dipole));
    return pre_fct * post_fct * au_to_barns * 1.0e-6;
}

double sigma_tot_spherical_symetry(const double& photon, const Eigen::Vector3cd& dipole) {
    double pre_fct  = 16. / 3. * std::pow(M_PI, 3) * photon / speed_of_light;
    double post_fct = dipole.squaredNorm();
    return pre_fct * post_fct * au_to_barns * 1.0e-6;
}


double photonEeV(const double k, const double ionization_pot) {
    double res = (k * k * 0.5 + ionization_pot) * 27.211385;
    return res;
}

double k_length(const double energy_eV, const double ionization_pot) {
    return sqrt(2.0 * (energy_eV / 27.211385 - ionization_pot));
}

Eigen::VectorXcd fetch_coulomb_wf(const int& lmax, const Eigen::Vector3d& kvec,
                                  const Eigen::VectorXd& norms) {
    std::vector<std::vector<std::vector<double>>> Dfact(lmax + 1);

    for (int l = 0; l <= lmax; l++) {
        Dfact[l].resize(l + 1);
        for (int p = 0; p <= l; p++)
            Dfact[l][p].resize(l - p + 1);
    }

    using std::pow;
    using std::sqrt;
    using boost::math::double_factorial;

    for (int l = 0; l <= lmax; l++)
        for (int p = 0; p <= l; p++)
            for (int q = 0; q <= l - p; q++) {
                Dfact[l][p][q] = 0.0;
                for (int m = -l; m <= l; m++)
                    Dfact[l][p][q] += Harmonics::NoNormCalcClmR(l, m, p, q, l - p - q) *
                                      Harmonics::Y(l, m, kvec);

                int r = l - p - q;
                int argx = (p == 0) ? 1 : 2 * p - 1;
                int argy = (q == 0) ? 1 : 2 * q - 1;
                int argz = (r == 0) ? 1 : 2 * r - 1;

                Dfact[l][p][q] *= pow(M_PI, 0.25) / pow(2.0, 0.25 + l) *
                                  sqrt(double_factorial<double>(argx) *
                                       double_factorial<double>(argy) *
                                       double_factorial<double>(argz));
            };

    std::vector<std::complex<double>> vec_cont =
            Gamess::order_set<std::complex<double>, double>(Dfact);

    Eigen::VectorXcd vec = Eigen::Map<Eigen::VectorXcd>(vec_cont.data(), vec_cont.size());
    for (int i = 0; i < vec.size(); i++)
        vec(i) /= norms.tail(vec.size())(i);

    return std::move(vec);
}
