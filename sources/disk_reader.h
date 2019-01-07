#ifndef DISK_READER_H
#define DISK_READER_H

#include <string>

#include <eigen3/Eigen/Dense>

#include "input_data.h"

class Disk_reader {
   public:
    Disk_reader(const Input_data &data);

    Eigen::MatrixXcd load_S(const std::string &path) const;
    Eigen::MatrixXcd load_H(const std::string &path) const;
    Eigen::MatrixXcd load_Dipx(const std::string &path) const;
    Eigen::MatrixXcd load_Dipy(const std::string &path) const;
    Eigen::MatrixXcd load_Dipz(const std::string &path) const;
    Eigen::MatrixXcd load_Gradx(const std::string &path) const;
    Eigen::MatrixXcd load_Grady(const std::string &path) const;
    Eigen::MatrixXcd load_Gradz(const std::string &path) const;

    Eigen::MatrixXd load_HFv(const std::string &path) const;
    Eigen::VectorXd load_HFe(const std::string &path) const;

    Eigen::VectorXd load_norms(const std::string &path) const;

   protected:
    Eigen::MatrixXcd load_matrix1E_bin(const std::string &path, const int &position) const;

   private:
    static constexpr int matrices1E_number = 20;

    int basis_l;
    int basis_lnk;
};

#endif