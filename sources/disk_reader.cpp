#include "disk_reader.h"

#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <complex>


Disk_reader::Disk_reader(const Input_data &data) {
    auto bnkl = std::stoi(data.first("NUMBER_GTO"));
    auto bkl = std::stoi(data.first("NUMBER_PWGTO"));
    
    basis_lnk = bnkl;    
    basis_l = bnkl + bkl;
}


Eigen::MatrixXcd Disk_reader::load_matrix1E_bin(const std::string &path, const int &position) const {
    std::ifstream file1E(path, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file1E.is_open())
        throw std::runtime_error("Unable to open 1E file.");

    auto bl_sqrt = basis_l * basis_l;

    std::streampos size1E = file1E.tellg();
    int complex_size      = size1E * sizeof(char) / sizeof(double);
    complex_size /= 2;

    if (complex_size != bl_sqrt * matrices1E_number)
        throw std::runtime_error("The size of 1E file does not match the basis. Have you used the correct 1E file?");

    std::vector<double> real(bl_sqrt);
    std::vector<double> imag(bl_sqrt);

    int chunk_size = size1E / matrices1E_number;
    file1E.seekg(position * chunk_size, std::ios::beg);
    chunk_size /= 2;
    file1E.read(reinterpret_cast<char *>(real.data()), chunk_size);
    file1E.read(reinterpret_cast<char *>(imag.data()), chunk_size);
    file1E.close();

    Eigen::MatrixXcd ints(basis_l, basis_l);

    std::transform(real.begin(), real.end(), imag.begin(), ints.data(),
                   [](double &dr, double &di) {
                       return std::complex<double>(dr, di);
                   });
    return ints.transpose();
}

Eigen::MatrixXcd Disk_reader::load_S(const std::string &path) const {
    return load_matrix1E_bin(path, 0);
}

Eigen::MatrixXcd Disk_reader::load_H(const std::string &path) const {
    return load_matrix1E_bin(path, 3);
}

Eigen::MatrixXcd Disk_reader::load_Dipx(const std::string &path) const {
    return load_matrix1E_bin(path, 4);
}

Eigen::MatrixXcd Disk_reader::load_Dipy(const std::string &path) const {
    return load_matrix1E_bin(path, 5);
}

Eigen::MatrixXcd Disk_reader::load_Dipz(const std::string &path) const {
    return load_matrix1E_bin(path, 6);
}

Eigen::MatrixXcd Disk_reader::load_Gradx(const std::string &path) const {
    return load_matrix1E_bin(path, 13);
}

Eigen::MatrixXcd Disk_reader::load_Grady(const std::string &path) const {
    return load_matrix1E_bin(path, 14);
}

Eigen::MatrixXcd Disk_reader::load_Gradz(const std::string &path) const {
    return load_matrix1E_bin(path, 15);
}

Eigen::MatrixXd Disk_reader::load_HFv(const std::string &path) const {

    Eigen::MatrixXd mat(basis_lnk, basis_lnk);

    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the HFv file.");

    auto size = basis_lnk * basis_lnk;
    for (int i = 0; i < size; ++i)
        file >> mat.data()[i];

    file.close();
    return mat;
}

Eigen::VectorXd Disk_reader::load_HFe(const std::string &path) const {
    Eigen::VectorXd vec(basis_lnk);

    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the HFv file.");

    for (int i = 0; i < basis_lnk; ++i)
        file >> vec(i);

    file.close();

    return vec;
}


Eigen::VectorXd Disk_reader::load_norms(const std::string &path) const {
    Eigen::VectorXd vec(basis_l);

    std::ifstream file(path, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the norms file.");

    std::streampos size = file.tellg();
    int double_size     = size * sizeof(char) / sizeof(double);
    if (double_size != basis_l)
        throw std::runtime_error("Size of the norms file is not consistent with basis. Have you used the correct norms file?");

    file.seekg(0, std::ios::beg);
    file.read(reinterpret_cast<char *>(vec.data()), size);
    file.close();

    return vec;
}
