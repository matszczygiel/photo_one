#ifndef FUN_H
#define FUN_H

#include <complex>

#include "constants.h"
#include "disk_reader.h"
#include "gamess.h"
#include "harmonics.h"

#include <eigen3/Eigen/Dense>

template <class type>
inline type scalar_prod(const Eigen::Matrix<type, Eigen::Dynamic, 1>& vec1,
                        const Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>& kernel,
                        const Eigen::Matrix<type, Eigen::Dynamic, 1>& vec2) {
    return vec1.dot(kernel * vec2);
}

double dsigma(const double& photon, const Eigen::Vector3d& polarization,const Eigen::Vector3cd& dipole );

double sigma_tot_spherical_symetry(const double& photon, const Eigen::Vector3cd& dipole);

double photonEeV(const double k, const double ionization_pot);

double k_length(const double energy_eV, const double ionization_pot);

Eigen::VectorXcd fetch_coulomb_wf(const int& lmax, const Eigen::Vector3d& kvec, const Eigen::VectorXd& norms);

#endif
