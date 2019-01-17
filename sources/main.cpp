/* A photoionization cros section evaluation program. The calculation is made by fixing asymptotic of the continuum wave function.
 * The final state is build as antisymetrized product of bound and continuum state. The mean field approximation is used to correct input orbitals.
 * The input is made by specyfing the asymptotic of continuum wave function.
 * Orthogonality of the continuum wave function to the ground state orbital is optional.
 * Two selection methods of orbitals are avaliable.
 *
 * The program uses the Eigen linear algebra libary.
 *
 * M. S. Szczygiel 2018
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "constants.h"
#include "disk_reader.h"
#include "functions.h"
#include "gamess.h"
#include "harmonics.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {
    auto start = chrono::system_clock::now();

    if (!(argc == 3 || argc == 2)) {
        cout << " Proper usage: ./photo <input name> <settings>\n";
        return EXIT_SUCCESS;
    }

    string input   = argv[1];
    string setting = "n";
    if (argc == 3)
        setting = argv[2];

    ifstream ifile(input);
    const Input_data data(ifile);
    ifile.close();

    const Disk_reader reader(data);

    const std::string enrg_i = data.first("PATH_IN") + data.first("FILE_HF_I_EN");
    const std::string orbs_i = data.first("PATH_IN") + data.first("FILE_HF_I_VEC");

    cout.precision(5);

    ////////////////////////////

    const auto energies_i = reader.load_HFe(enrg_i);
    const auto orbitals_i = reader.load_HFv(orbs_i);

    const double photon   = std::stof(data.first("PHOTON_EN")) / au_to_ev;
    const double en_final = energies_i[0] + photon;
    if (en_final < 0)
        throw std::runtime_error("Energy below the ionization threshold!");

    const double kval = sqrt(2. * en_final);

    ///////////////////////////////

    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << kval;
    const std::string k_str = stream.str();
    std::cout << "K:  " << k_str << "\n";

    if (setting == "-dump")
        return 0;

    const auto lmax     = std::stoi(data.first("MAX_L"));
    const double ptheta = std::stof(data.first("POL_THETA"));
    const double pphi   = std::stof(data.first("POL_PHI"));
    /////////////
    vector<string> norm_files, ints_files, Rints_files;
    vector<double> theta;

    const int job_size = data.size("FILES_NORM") - 1;

    for (int i = 0; i < job_size; ++i) {
        norm_files.push_back(data.first("PATH_IN") + data.first("FILES_NORM") + k_str + data("FILES_NORM", i + 1));
        ints_files.push_back(data.first("PATH_IN") + data.first("FILES_1E") + k_str + data("FILES_1E", i + 1));
        theta.push_back(std::stof(data("K_THETA", i)));
    }

    double phi = std::stof(data.first("K_PHI"));

    vector<double> sigma(job_size);
    std::string gauge = data.first("GAUGE");
    std::transform(gauge.begin(), gauge.end(), gauge.begin(), ::tolower);

    const auto bnkl = std::stoi(data.first("NUMBER_GTO"));
    const auto bkl  = std::stoi(data.first("NUMBER_PWGTO"));

    for (int i = 0; i < job_size; ++i) {
        Vector3d kvec;
        kvec(0) = kval * sin(theta.at(i)) * cos(phi);
        kvec(1) = kval * sin(theta.at(i)) * sin(phi);
        kvec(2) = kval * cos(theta.at(i));

        const auto norms = reader.load_norms(norm_files.at(i));

        const auto vecc = fetch_coulomb_wf(lmax, kvec, norms);
        const auto veci = orbitals_i.col(0);
        //////////////////////

        MatrixXcd Dx, Dy, Dz;
        if (gauge == "dipole") {
            Dx = reader.load_Dipx(ints_files.at(i));
            Dy = reader.load_Dipy(ints_files.at(i));
            Dz = reader.load_Dipz(ints_files.at(i));
        } else if (gauge == "velocity") {
            Dx = reader.load_Gradx(ints_files.at(i)) / photon;
            Dy = reader.load_Grady(ints_files.at(i)) / photon;
            Dz = reader.load_Gradz(ints_files.at(i)) / photon;
        } else
            throw std::runtime_error("Invalid argument for GAUGE.");

        const auto S = reader.load_S(ints_files.at(i));

        // Dipole moment
        Vector3cd T;

        T(0) = vecc.dot(Dx.bottomLeftCorner(bkl, bnkl) * veci);
        T(1) = vecc.dot(Dy.bottomLeftCorner(bkl, bnkl) * veci);
        T(2) = vecc.dot(Dz.bottomLeftCorner(bkl, bnkl) * veci);

        //normalize to energy scale
        T *= std::sqrt(kval);

        cout << " Dipole moment: \n";
        cout << T << "\n\n";

        Vector3d j;
        j(0) = sin(ptheta) * cos(pphi);
        j(1) = sin(ptheta) * sin(pphi);
        j(2) = cos(ptheta);

        // keep this for H
        sigma.at(i) = sigma_tot_spherical_symetry(photon, T);

        //keep this for H2+
        //sigma.at(i).at(k) += dsigma(photon, j, T);

        cout << " Photon energy [eV]:  " << fixed << setprecision(3) << data.first("PHOTON_EN") << "\n";
        cout << " Cross section :      " << fixed << setprecision(4) << sigma.at(i) << "\n";
        cout << " \n\n\n";
    }

    //write results

    bool write;
    char token;
    std::string arg = "WRITE";

    token = std::tolower(*(data.first(arg).begin()));
    if (token == 'y')
        write = true;
    else if (token == 'n')
        write = false;
    else
        throw std::runtime_error("Invalid argument for" + arg + ".");

    if (write) {
        string res_path = data.first("PATH_OUT") + data.first("FILE_OUT");

        std::ofstream outfile(res_path, std::ios_base::app);
        outfile << std::fixed;
        outfile << "****** " << data.first("NAME") << " ******\n";
        outfile << "Photon [eV]         " << data.first("PHOTON_EN") << "\n";
        outfile << "k ( phi)            " << std::setprecision(3) << phi << "\n";
        outfile << "j (theta, phi)      " << std::setprecision(3) << ptheta << "\t" << std::setprecision(3) << pphi << "\n";
        outfile << "========\n";
        outfile << "k(theta)\ttot sigma\n";
        for (int i = 0; i < job_size; ++i) {
            outfile << std::setprecision(3) << theta.at(i) << "\t\t";
            outfile << std::setprecision(5) << sigma.at(i) << "\n";
        }
        outfile << "\n";
        outfile.close();
    }

    // successful exit
    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;
    cout << " CPU time: " << setprecision(5) << fixed << elapsed_seconds.count() << " s\n";
    cout << "=========================================================================="
         << "\n"
         << "\n";

    return EXIT_SUCCESS;
}
