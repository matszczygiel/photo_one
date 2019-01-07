#pragma once

#include <eigen3/Eigen/Dense>

namespace Harmonics
{
/**
         * Returns the coefficents of cartesian monomial in r^l * Y_{lm}, where Y_{lm} is REAL spherical harmonic.
         */
double NoNormCalcClmR(const int &l, const int &m, const int &lx, const int &ly, const int &lz);

/**
         * Return value of r^l * Y_{lm} (r), where Y_{lm} is REAL spherical harmonic.
         */
double NoNormYrl(const int &l, const int &m, const Eigen::Vector3d &r);

/**
         * Return value of Y_{lm} (r), where Y_{lm} is REAL spherical harmonic.
         */
double Y(const int &l, const int &m, const Eigen::Vector3d &r);
};

