//
//  Dxx_coeffs.h
//  magpie
//
//  Created by admin on 08/09/2025.
//

#ifndef Dxx_coeffs_h
#define Dxx_coeffs_h

#include <vector>
#include <cmath>

std::vector<double> D00_coeffs_x(double K0y, double R0y, double Kx0, double Rx0, double h, double D,  double nu){
    /*param K0y:
    :param R0y:
    :param Kx0:
    :param Rx0:
    :param h:
    :param D:
    :param nu:
    :return*/
    double u00_c = ((2 * K0y * Rx0 * (std::pow(R0y, 2)) * (std::pow(h, 6)) + 4 * K0y * D * (std::pow(R0y, 2)) * (std::pow(h, 5)) + 8 * K0y * Rx0 * D * R0y * (
            std::pow(h, 5)) - 8 * K0y * (std::pow(D, 2)) * R0y * (std::pow(h, 4)) * (std::pow(nu, 2)) + 16 * K0y * (std::pow(D, 2)) * R0y * (
                      std::pow(h, 4)) - 12 * Rx0 * (std::pow(D, 2)) * R0y * (std::pow(h, 2)) * (std::pow(nu, 2)) + 24 * Rx0 * (std::pow(D, 2)) * R0y * (
                      std::pow(h, 2)) * nu + 24 * Rx0 * (std::pow(D, 2)) * R0y * (std::pow(h, 2)) + 16 * (std::pow(D, 3)) * R0y * h * (
                      std::pow(nu, 3)) - 56 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) + 48 * (std::pow(D, 3)) * R0y * h + 8 * K0y * Rx0 * (
                      std::pow(D, 2)) * (std::pow(h, 4)) - 16 * K0y * (std::pow(D, 3)) * (std::pow(h, 3)) * (std::pow(nu, 2)) + 16 * K0y * (std::pow(D, 3)) * (
                      std::pow(h, 3)) - 24 * Rx0 * (std::pow(D, 3)) * h * (std::pow(nu, 2)) + 48 * Rx0 * (std::pow(D, 3)) * h * nu + 48 * Rx0 * (
                      std::pow(D, 3)) * h + 16 * (std::pow(D, 4)) * (std::pow(nu, 4)) - 112 * (std::pow(D, 4)) * (std::pow(nu, 2)) + 96 * (std::pow(D, 4))) / (
                     D * (2 * D + R0y * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - (
                     8 * (- 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * Rx0 * h * D * nu + 8 * (std::pow(D, 2)) + 4 * Rx0 * h * D)) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) - (4 * D * nu) / (2 * D + R0y * h) - (4 * D * nu) / (2 * D + Rx0 * h) - (
                     2 * (8 * (std::pow(D, 2)) * nu + 2 * D * R0y * h * nu + 2 * D * Rx0 * h * nu)) / (
                     (2 * D + R0y * h) * (2 * D + Rx0 * h)) - (
                     8 * (- 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * R0y * h * D * nu + 8 * (std::pow(D, 2)) + 4 * R0y * h * D)) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + (2 * Kx0 * R0y * (std::pow(Rx0, 2)) * (std::pow(h, 6)) + 4 * Kx0 * D * (std::pow(Rx0, 2)) * (
            std::pow(h, 5)) + 8 * Kx0 * R0y * D * Rx0 * (std::pow(h, 5)) - 8 * Kx0 * (std::pow(D, 2)) * Rx0 * (std::pow(h, 4)) * (
                                         std::pow(nu, 2)) + 16 * Kx0 * (std::pow(D, 2)) * Rx0 * (std::pow(h, 4)) - 12 * R0y * (
                                         std::pow(D, 2)) * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 24 * R0y * (
                                         std::pow(D, 2)) * Rx0 * (std::pow(h, 2)) * nu + 24 * R0y * (std::pow(D, 2)) * Rx0 * (
                                         std::pow(h, 2)) + 16 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 3)) - 56 * (
                                         std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) + 48 * (
                                         std::pow(D, 3)) * Rx0 * h + 8 * Kx0 * R0y * (std::pow(D, 2)) * (
                                         std::pow(h, 4)) - 16 * Kx0 * (std::pow(D, 3)) * (std::pow(h, 3)) * (std::pow(nu, 2)) + 16 * Kx0 * (
                                         std::pow(D, 3)) * (std::pow(h, 3)) - 24 * R0y * (std::pow(D, 3)) * h * (
                                         std::pow(nu, 2)) + 48 * R0y * (std::pow(D, 3)) * h * nu + 48 * R0y * (
                                         std::pow(D, 3)) * h + 16 * (std::pow(D, 4)) * (std::pow(nu, 4)) - 112 * (std::pow(D, 4)) * (
                                         std::pow(nu, 2)) + 96 * (std::pow(D, 4))) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) + 20);
    double u10_c = ((2 * (4 * D + 4 * D * nu)) / (2 * D + Rx0 * h) - (8 * (
            4 * (std::pow(D, 2)) * (std::pow(nu, 2)) - 4 * (std::pow(D, 2)) + 2 * D * R0y * h - 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + (
                     2 * (8 * (std::pow(D, 2)) * nu + 8 * (std::pow(D, 2)) + 4 * D * R0y * h + 4 * D * R0y * h * nu)) / (
                     (2 * D + R0y * h) * (2 * D + Rx0 * h)) + (
                     32 * (std::pow(D, 4)) * nu - 96 * (std::pow(D, 4)) + 96 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) - 48 * (std::pow(D, 3)) * R0y * h - 48 * (std::pow(D, 3)) * Rx0 * h + 16 * (
                             std::pow(D, 3)) * R0y * h * nu + 16 * (std::pow(D, 3)) * Rx0 * h * nu - 24 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 48 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 16 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 3)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (
                     D * (2 * D + R0y * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - (
                     32 * (std::pow(D, 4)) * nu + 64 * (std::pow(D, 4)) - 96 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) + 32 * (std::pow(D, 4)) * (std::pow(nu, 4)) + 32 * (std::pow(D, 3)) * R0y * h + 32 * (
                             std::pow(D, 3)) * Rx0 * h + 64 * (std::pow(D, 3)) * R0y * h * nu + 16 * (
                             std::pow(D, 3)) * Rx0 * h * nu + 16 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) - 32 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 16 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 16 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 32 * (std::pow(D, 2)) * R0y * Rx0 * (
                             std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) + (
                     32 * D * R0y * h * nu) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) - 8);
    double u20_c = ((Rx0 * D * (std::pow(R0y, 2)) * (std::pow(h, 3)) + 2 * (std::pow(D, 2)) * (std::pow(R0y, 2)) * (std::pow(h, 2)) + 4 * Rx0 * (std::pow(D, 2)) * R0y * (
            std::pow(h, 2)) - 4 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) + 8 * (std::pow(D, 3)) * R0y * h + 4 * Rx0 * (std::pow(D, 3)) * h - 8 * (
                      std::pow(D, 4)) * (std::pow(nu, 2)) + 8 * (std::pow(D, 4))) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     2 * (4 * nu * (std::pow(D, 2)) + 2 * R0y * h * nu * D)) / ((2 * D + R0y * h) * (2 * D + Rx0 * h)) - (
                     4 * D * nu) / (2 * D + Rx0 * h) + (
                     32 * (std::pow(D, 4)) * nu - 16 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (std::pow(nu, 3)) + 16 * (std::pow(D, 4)) * (
                     std::pow(nu, 4)) + 16 * (std::pow(D, 3)) * R0y * h * nu + 16 * (std::pow(D, 3)) * Rx0 * h * nu - 8 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 8 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 4 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                             std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) + 1);
    double u01_c = ((2 * (4 * D + 4 * D * nu)) / (2 * D + R0y * h) - (8 * (
            4 * (std::pow(D, 2)) * (std::pow(nu, 2)) - 4 * (std::pow(D, 2)) - 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + (
                     2 * (8 * (std::pow(D, 2)) * nu + 8 * (std::pow(D, 2)) + 4 * D * Rx0 * h + 4 * D * Rx0 * h * nu)) / (
                     (2 * D + R0y * h) * (2 * D + Rx0 * h)) + (
                     32 * (std::pow(D, 4)) * nu - 96 * (std::pow(D, 4)) + 96 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) - 48 * (std::pow(D, 3)) * R0y * h - 48 * (std::pow(D, 3)) * Rx0 * h + 16 * (
                             std::pow(D, 3)) * R0y * h * nu + 16 * (std::pow(D, 3)) * Rx0 * h * nu - 24 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 48 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 16 * (
                             std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 3)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (
                     D * (2 * D + Rx0 * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - (
                     32 * (std::pow(D, 4)) * nu + 64 * (std::pow(D, 4)) - 96 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) + 32 * (std::pow(D, 4)) * (std::pow(nu, 4)) + 32 * (std::pow(D, 3)) * R0y * h + 32 * (
                             std::pow(D, 3)) * Rx0 * h + 16 * (std::pow(D, 3)) * R0y * h * nu + 64 * (
                             std::pow(D, 3)) * Rx0 * h * nu + 16 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) - 16 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 32 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 16 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 32 * (std::pow(D, 2)) * R0y * Rx0 * (
                             std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) + (
                     32 * D * Rx0 * h * nu) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) - 8);
    double u02_c = ((R0y * D * (std::pow(Rx0, 2)) * (std::pow(h, 3)) + 2 * (std::pow(D, 2)) * (std::pow(Rx0, 2)) * (std::pow(h, 2)) + 4 * R0y * (std::pow(D, 2)) * Rx0 * (
            std::pow(h, 2)) - 4 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) + 8 * (std::pow(D, 3)) * Rx0 * h + 4 * R0y * (std::pow(D, 3)) * h - 8 * (
                      std::pow(D, 4)) * (std::pow(nu, 2)) + 8 * (std::pow(D, 4))) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     2 * (4 * nu * (std::pow(D, 2)) + 2 * Rx0 * h * nu * D)) / ((2 * D + R0y * h) * (2 * D + Rx0 * h)) - (
                     4 * D * nu) / (2 * D + R0y * h) + (
                     32 * (std::pow(D, 4)) * nu - 16 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (std::pow(nu, 3)) + 16 * (std::pow(D, 4)) * (
                     std::pow(nu, 4)) + 16 * (std::pow(D, 3)) * R0y * h * nu + 16 * (std::pow(D, 3)) * Rx0 * h * nu - 8 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 8 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 4 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                             std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) + 1);
    double u11_c = (2 - (2 * (2 * D - Rx0 * h)) / (2 * D + Rx0 * h) - (
            2 * (12 * (std::pow(D, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h - R0y * Rx0 * (std::pow(h, 2)))) / (
                     (2 * D + R0y * h) * (2 * D + Rx0 * h)) - (
                     32 * (std::pow(D, 4)) * nu - 64 * (std::pow(D, 4)) + 64 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) - 32 * (std::pow(D, 3)) * R0y * h - 32 * (std::pow(D, 3)) * Rx0 * h + 16 * (
                             std::pow(D, 3)) * R0y * h * nu + 16 * (std::pow(D, 3)) * Rx0 * h * nu - 16 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (
                     D * (2 * D + R0y * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - (
                     32 * (std::pow(D, 4)) * nu - 64 * (std::pow(D, 4)) + 64 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) - 32 * (std::pow(D, 3)) * R0y * h - 32 * (std::pow(D, 3)) * Rx0 * h + 16 * (
                             std::pow(D, 3)) * R0y * h * nu + 16 * (std::pow(D, 3)) * Rx0 * h * nu - 16 * (
                             std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (
                     D * (2 * D + Rx0 * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - (2 * (2 * D - R0y * h)) / (2 * D + R0y * h));

    return {u00_c, u10_c, u20_c, u01_c, u02_c, u11_c};
}

std::vector<double> D01_coeffs_x(double K0y, double R0y, double Rx0, double h, double D, double nu){
    /*param K0y:
    :param R0y:
    :param Rx0:
    :param h:
    :param D:
    :param nu:
    :return*/
    double u01_c = ((4 * std::pow(D, 2) * std::pow(nu, 2) - 4 * std::pow(D, 2) - 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * std::pow(h, 2)) / (4 * std::pow(D, 2) - 4 * std::pow(D, 2) * std::pow(nu, 2) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * std::pow(h, 2)) - (8 * (4 * D + 4 * D * nu)) / (2 * D + R0y * h) - (4 * D * nu) / (2 * D + R0y * h) + (2 * K0y * Rx0 * std::pow(R0y, 2) * std::pow(h, 6) + 4 * K0y * D * std::pow(R0y, 2) * std::pow(h, 5) + 8 * K0y * Rx0 * D * R0y * std::pow(h, 5) - 8 * K0y * std::pow(D, 2) * R0y * std::pow(h, 4) * std::pow(nu, 2) + 16 * K0y * std::pow(D, 2) * R0y * std::pow(h, 4) - 14 * Rx0 * std::pow(D, 2) * R0y * std::pow(h, 2) * std::pow(nu, 2) + 28 * Rx0 * std::pow(D, 2) * R0y * std::pow(h, 2) * nu + 24 * Rx0 * std::pow(D, 2) * R0y * std::pow(h, 2) - 20 * std::pow(D, 3) * R0y * h * std::pow(nu, 2) + 40 * std::pow(D, 3) * R0y * h * nu + 48 * std::pow(D, 3) * R0y * h + 8 * K0y * Rx0 * std::pow(D, 2) * std::pow(h, 4) - 16 * K0y * std::pow(D, 3) * std::pow(h, 3) * std::pow(nu, 2) + 16 * K0y * std::pow(D, 3) * std::pow(h, 3) - 28 * Rx0 * std::pow(D, 3) * h * std::pow(nu, 2) + 56 * Rx0 * std::pow(D, 3) * h * nu + 48 * Rx0 * std::pow(D, 3) * h + 40 * std::pow(D, 4) * std::pow(nu, 4) - 80 * std::pow(D, 4) * std::pow(nu, 3) - 136 * std::pow(D, 4) * std::pow(nu, 2) + 80 * std::pow(D, 4) * nu + 96 * std::pow(D, 4)) / (D * (2 * D + R0y * h) * (4 * std::pow(D, 2) - 4 * std::pow(D, 2) * std::pow(nu, 2) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * std::pow(h, 2))) - (8 * D * Rx0 * h * nu) / (4 * std::pow(D, 2) - 4 * std::pow(D, 2) * std::pow(nu, 2) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * std::pow(h, 2)) + 20) ;
//    double u01_c = ((4 * (std::pow(D, 2)) * (std::pow(nu, 2)) - 4 * (std::pow(D, 2)) - 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2))) / (
//            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2))) - (
//                     8 * (4 * D + 4 * D * nu)) / (2 * D + R0y * h) - (4 * D * nu) / (2 * D + R0y * h) + (
//                     2 * K0y * Rx0 * (std::pow(R0y, 2)) * (std::pow(h, 6)) + 4 * K0y * D * (std::pow(R0y, 2)) * (
//                     std::pow(h, 5)) + 8 * K0y * Rx0 * D * R0y * (std::pow(h, 5)) - 8 * K0y * (std::pow(D, 2)) * R0y * (std::pow(h, 4)) * (
//                             std::pow(nu, 2)) + 16 * K0y * (std::pow(D, 2)) * R0y * (std::pow(h, 4)) - 14 * Rx0 * (std::pow(D, 2)) * R0y * (
//                             std::pow(h, 2)) * (std::pow(nu, 2)) + 28 * Rx0 * (std::pow(D, 2)) * R0y * (std::pow(h, 2)) * nu + 24 * Rx0 * (
//                             std::pow(D, 2)) * R0y * (std::pow(h, 2)) - 20 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) + 40 * (
//                             std::pow(D, 3)) * R0y * h * nu + 48 * (std::pow(D, 3)) * R0y * h + 8 * K0y * Rx0 * (std::pow(D, 2)) * (
//                             std::pow(h, 4)) - 16 * K0y * (std::pow(D, 3)) * (std::pow(h, 3)) * (std::pow(nu, 2)) + 16 * K0y * (std::pow(D, 3)) * (
//                             std::pow(h, 3)) - 28 * Rx0 * (std::pow(D, 3)) * h * (std::pow(nu, 2)) + 56 * Rx0 * (
//                             std::pow(D, 3)) * h * nu + 48 * Rx0 * (std::pow(D, 3)) * h + 40 * (std::pow(D, 4)) * (std::pow(nu, 4)) - 80 * (
//                             std::pow(D, 4)) * (std::pow(nu, 3)) - 136 * (std::pow(D, 4)) * (std::pow(nu, 2)) + 80 * (std::pow(D, 4)) * nu + 96 * (
//                             std::pow(D, 4))) / (D * (2 * D + R0y * h) * (
//            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
//                     8 * D * Rx0 * h * nu) / (
//                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
//                     std::pow(h, 2))) + 20);

    double u11_c = ((8 * (2 * D - R0y * h)) / (2 * D + R0y * h) + (
            32 * (std::pow(D, 4)) * nu - 96 * (std::pow(D, 4)) + 96 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (std::pow(nu, 3)) - 48 * (
            std::pow(D, 3)) * R0y * h - 48 * (std::pow(D, 3)) * Rx0 * h + 16 * (std::pow(D, 3)) * R0y * h * nu + 16 * (
                    std::pow(D, 3)) * Rx0 * h * nu - 24 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                    std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) - 8);

    double u21_c = ((Rx0 * D * (std::pow(R0y, 2)) * (std::pow(h, 3)) + 2 * (std::pow(D, 2)) * (std::pow(R0y, 2)) * (std::pow(h, 2)) + 4 * Rx0 * (std::pow(D, 2)) * R0y * (
            std::pow(h, 2)) - 4 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) + 8 * (std::pow(D, 3)) * R0y * h + 4 * Rx0 * (std::pow(D, 3)) * h - 8 * (
                      std::pow(D, 4)) * (std::pow(nu, 2)) + 8 * (std::pow(D, 4))) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) + 1);

    double u00_c = ((- 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * R0y * h * D * nu + 8 * (std::pow(D, 2)) + 4 * R0y * h * D) / (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2))) + (
                     2 * (- 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * Rx0 * h * D * nu + 8 * (std::pow(D, 2)) + 4 * Rx0 * h * D)) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + (16 * D * nu) / (2 * D + R0y * h) - (
                     32 * (std::pow(D, 4)) * nu + 32 * (std::pow(D, 4)) - 48 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) + 16 * (std::pow(D, 4)) * (std::pow(nu, 4)) + 16 * (std::pow(D, 3)) * R0y * h + 16 * (
                             std::pow(D, 3)) * Rx0 * h + 16 * (std::pow(D, 3)) * R0y * h * nu + 32 * (
                             std::pow(D, 3)) * Rx0 * h * nu + 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) - 24 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) + 8 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 3)) - 16 * (
                             std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (
                             std::pow(nu, 2)) + 16 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) - 8);

    double u02_c = ((2 * (4 * D + 4 * D * nu)) / (2 * D + R0y * h) + (16 * D * nu) / (2 * D + R0y * h) - (
            64 * (std::pow(D, 4)) * nu + 32 * (std::pow(D, 4)) - 64 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 64 * (std::pow(D, 4)) * (std::pow(nu, 3)) + 32 * (
            std::pow(D, 4)) * (std::pow(nu, 4)) + 16 * (std::pow(D, 3)) * R0y * h + 16 * (std::pow(D, 3)) * Rx0 * h + 32 * (
                    std::pow(D, 3)) * R0y * h * nu + 32 * (std::pow(D, 3)) * Rx0 * h * nu + 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                    std::pow(h, 2)) - 16 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 16 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 8 * (
                    std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 16 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (
                     D * (2 * D + R0y * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - 8);

    double u03_c = ((16 * (std::pow(D, 4)) * nu - 8 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 16 * (std::pow(D, 4)) * (std::pow(nu, 3)) + 8 * (std::pow(D, 4)) * (
            std::pow(nu, 4)) + 8 * (std::pow(D, 3)) * R0y * h * nu + 8 * (std::pow(D, 3)) * Rx0 * h * nu - 4 * (std::pow(D, 3)) * R0y * h * (
                      std::pow(nu, 2)) - 4 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 2 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (
                      std::pow(nu, 2)) + 4 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     4 * D * nu) / (2 * D + R0y * h) + 1);

    double u12_c = (2 - (16 * (std::pow(D, 4)) * nu - 32 * (std::pow(D, 4)) + 32 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 16 * (std::pow(D, 4)) * (std::pow(nu, 3)) - 16 * (
            std::pow(D, 3)) * R0y * h - 16 * (std::pow(D, 3)) * Rx0 * h + 8 * (std::pow(D, 3)) * R0y * h * nu + 8 * (
                          std::pow(D, 3)) * Rx0 * h * nu - 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 4 * (
                          std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     2 * (2 * D - R0y * h)) / (2 * D + R0y * h));

    double u10_c = ((2 * (
            4 * (std::pow(D, 2)) * (std::pow(nu, 2)) - 4 * (std::pow(D, 2)) + 2 * D * R0y * h - 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) - (
                     16 * (std::pow(D, 4)) * nu - 32 * (std::pow(D, 4)) + 32 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 16 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) - 16 * (std::pow(D, 3)) * R0y * h - 16 * (std::pow(D, 3)) * Rx0 * h + 8 * (
                             std::pow(D, 3)) * R0y * h * nu + 8 * (std::pow(D, 3)) * Rx0 * h * nu - 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                             std::pow(h, 2)) + 16 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 8 * (std::pow(D, 3)) * R0y * h * (
                             std::pow(nu, 3)) + 4 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + R0y * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     4 * D * R0y * h * nu) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + 2);

    return {u01_c, u11_c, u21_c, u00_c, u02_c, u03_c, u12_c, u10_c};
}
std::vector<double> D02_coeffs_x(double K0y, double R0y, double h, double D, double nu){
    /*param K0y:
    :param R0y:
    :param h:
    :param D:
    :param nu:
    :return*/
    double u02_c = ((2 * K0y * R0y * (std::pow(h, 4)) + 4 * K0y * D * (std::pow(h, 3)) - 12 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 24 * (
            std::pow(D, 2)) * nu + 24 * (std::pow(D, 2))) / (D * (2 * D + R0y * h)) - (8 * D * nu) / (2 * D + R0y * h) - (
                     8 * (4 * D + 4 * D * nu)) / (2 * D + R0y * h) + 20);

    double u12_c = ((8 * (2 * D - R0y * h)) / (2 * D + R0y * h) + (8 * (std::pow(D, 2)) * nu - 24 * (std::pow(D, 2))) / (
            D * (2 * D + R0y * h)) - 8);

    double u22_c = ((2 * (std::pow(D, 2)) + R0y * h * D) / (D * (2 * D + R0y * h)) + 1);

    double u01_c = ((2 * (4 * D + 4 * D * nu)) / (2 * D + R0y * h) - (
            - 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 16 * (std::pow(D, 2)) * nu + 8 * (std::pow(D, 2))) / (D * (2 * D + R0y * h)) + (
                     16 * D * nu) / (2 * D + R0y * h) - 8);
    double u03_c = u01_c;

    double u04_c = ((- 2 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * (std::pow(D, 2)) * nu) / (D * (2 * D + R0y * h)) - (4 * D * nu) / (
            2 * D + R0y * h) + 1);
    double u00_c = u04_c;

    double u13_c = (2 - (4 * (std::pow(D, 2)) * nu - 8 * (std::pow(D, 2))) / (D * (2 * D + R0y * h)) - (2 * (2 * D - R0y * h)) / (
            2 * D + R0y * h));
    double u11_c = u13_c;

    return {u02_c, u12_c, u22_c, u01_c, u03_c, u04_c, u00_c, u13_c, u11_c};
}

std::vector<double> D10_coeffs_x(double R0y, double Kx0, double Rx0, double h, double D, double nu){
    /*param R0y:
    :param Kx0:
    :param Rx0:
    :param h:
    :param D:
    :param nu:
    :return*/
    double u10_c = ((4 * (std::pow(D, 2)) * (std::pow(nu, 2)) - 4 * (std::pow(D, 2)) + 2 * D * R0y * h - 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2))) / (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2))) - (
                     8 * (4 * D + 4 * D * nu)) / (2 * D + Rx0 * h) - (4 * D * nu) / (2 * D + Rx0 * h) + (
                     2 * Kx0 * R0y * (std::pow(Rx0, 2)) * (std::pow(h, 6)) + 4 * Kx0 * D * (std::pow(Rx0, 2)) * (
                     std::pow(h, 5)) + 8 * Kx0 * R0y * D * Rx0 * (std::pow(h, 5)) - 8 * Kx0 * (std::pow(D, 2)) * Rx0 * (std::pow(h, 4)) * (
                             std::pow(nu, 2)) + 16 * Kx0 * (std::pow(D, 2)) * Rx0 * (std::pow(h, 4)) - 14 * R0y * (std::pow(D, 2)) * Rx0 * (
                             std::pow(h, 2)) * (std::pow(nu, 2)) + 28 * R0y * (std::pow(D, 2)) * Rx0 * (std::pow(h, 2)) * nu + 24 * R0y * (
                             std::pow(D, 2)) * Rx0 * (std::pow(h, 2)) - 20 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) + 40 * (
                             std::pow(D, 3)) * Rx0 * h * nu + 48 * (std::pow(D, 3)) * Rx0 * h + 8 * Kx0 * R0y * (std::pow(D, 2)) * (
                             std::pow(h, 4)) - 16 * Kx0 * (std::pow(D, 3)) * (std::pow(h, 3)) * (std::pow(nu, 2)) + 16 * Kx0 * (std::pow(D, 3)) * (
                             std::pow(h, 3)) - 28 * R0y * (std::pow(D, 3)) * h * (std::pow(nu, 2)) + 56 * R0y * (
                             std::pow(D, 3)) * h * nu + 48 * R0y * (std::pow(D, 3)) * h + 40 * (std::pow(D, 4)) * (std::pow(nu, 4)) - 80 * (
                             std::pow(D, 4)) * (std::pow(nu, 3)) - 136 * (std::pow(D, 4)) * (std::pow(nu, 2)) + 80 * (std::pow(D, 4)) * nu + 96 * (
                             std::pow(D, 4))) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     8 * D * R0y * h * nu) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + 20);

    double u20_c = ((2 * (4 * D + 4 * D * nu)) / (2 * D + Rx0 * h) + (16 * D * nu) / (2 * D + Rx0 * h) - (
            64 * (std::pow(D, 4)) * nu + 32 * (std::pow(D, 4)) - 64 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 64 * (std::pow(D, 4)) * (std::pow(nu, 3)) + 32 * (
            std::pow(D, 4)) * (std::pow(nu, 4)) + 16 * (std::pow(D, 3)) * R0y * h + 16 * (std::pow(D, 3)) * Rx0 * h + 32 * (
                    std::pow(D, 3)) * R0y * h * nu + 32 * (std::pow(D, 3)) * Rx0 * h * nu + 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                    std::pow(h, 2)) - 16 * (std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 16 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 8 * (
                    std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (std::pow(nu, 2)) + 16 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (
                     D * (2 * D + Rx0 * h) * (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2)))) - 8);

    double u30_c = ((16 * (std::pow(D, 4)) * nu - 8 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 16 * (std::pow(D, 4)) * (std::pow(nu, 3)) + 8 * (std::pow(D, 4)) * (
            std::pow(nu, 4)) + 8 * (std::pow(D, 3)) * R0y * h * nu + 8 * (std::pow(D, 3)) * Rx0 * h * nu - 4 * (std::pow(D, 3)) * R0y * h * (
                      std::pow(nu, 2)) - 4 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 2 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (
                      std::pow(nu, 2)) + 4 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     4 * D * nu) / (2 * D + Rx0 * h) + 1);

    double u00_c = ((2 * (- 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * R0y * h * D * nu + 8 * (std::pow(D, 2)) + 4 * R0y * h * D)) / (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2))) + (
                     - 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * Rx0 * h * D * nu + 8 * (std::pow(D, 2)) + 4 * Rx0 * h * D) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + (16 * D * nu) / (2 * D + Rx0 * h) - (
                     32 * (std::pow(D, 4)) * nu + 32 * (std::pow(D, 4)) - 48 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) + 16 * (std::pow(D, 4)) * (std::pow(nu, 4)) + 16 * (std::pow(D, 3)) * R0y * h + 16 * (
                             std::pow(D, 3)) * Rx0 * h + 32 * (std::pow(D, 3)) * R0y * h * nu + 16 * (
                             std::pow(D, 3)) * Rx0 * h * nu + 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) - 16 * (
                             std::pow(D, 3)) * R0y * h * (std::pow(nu, 2)) - 24 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) + 8 * (
                             std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 3)) - 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * (
                             std::pow(nu, 2)) + 16 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) - 8);

    double u11_c = ((8 * (2 * D - Rx0 * h)) / (2 * D + Rx0 * h) + (
            32 * (std::pow(D, 4)) * nu - 96 * (std::pow(D, 4)) + 96 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 32 * (std::pow(D, 4)) * (std::pow(nu, 3)) - 48 * (
            std::pow(D, 3)) * R0y * h - 48 * (std::pow(D, 3)) * Rx0 * h + 16 * (std::pow(D, 3)) * R0y * h * nu + 16 * (
                    std::pow(D, 3)) * Rx0 * h * nu - 24 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                    std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) - 8);

    double u12_c = ((R0y * D * (std::pow(Rx0, 2)) * (std::pow(h, 3)) + 2 * (std::pow(D, 2)) * (std::pow(Rx0, 2)) * (std::pow(h, 2)) + 4 * R0y * (std::pow(D, 2)) * Rx0 * (
            std::pow(h, 2)) - 4 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) + 8 * (std::pow(D, 3)) * Rx0 * h + 4 * R0y * (std::pow(D, 3)) * h - 8 * (
                      std::pow(D, 4)) * (std::pow(nu, 2)) + 8 * (std::pow(D, 4))) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
            std::pow(h, 2)))) + 1);

    double u21_c = (2 - (16 * (std::pow(D, 4)) * nu - 32 * (std::pow(D, 4)) + 32 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 16 * (std::pow(D, 4)) * (std::pow(nu, 3)) - 16 * (
            std::pow(D, 3)) * R0y * h - 16 * (std::pow(D, 3)) * Rx0 * h + 8 * (std::pow(D, 3)) * R0y * h * nu + 8 * (
                          std::pow(D, 3)) * Rx0 * h * nu - 8 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) + 4 * (
                          std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     2 * (2 * D - Rx0 * h)) / (2 * D + Rx0 * h));

    double u01_c = ((2 * (
            4 * (std::pow(D, 2)) * (std::pow(nu, 2)) - 4 * (std::pow(D, 2)) - 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) - (
                     16 * (std::pow(D, 4)) * nu - 32 * (std::pow(D, 4)) + 32 * (std::pow(D, 4)) * (std::pow(nu, 2)) - 16 * (std::pow(D, 4)) * (
                     std::pow(nu, 3)) - 16 * (std::pow(D, 3)) * R0y * h - 16 * (std::pow(D, 3)) * Rx0 * h + 8 * (
                             std::pow(D, 3)) * R0y * h * nu + 8 * (std::pow(D, 3)) * Rx0 * h * nu - 8 * (std::pow(D, 2)) * R0y * Rx0 * (
                             std::pow(h, 2)) + 16 * (std::pow(D, 3)) * Rx0 * h * (std::pow(nu, 2)) - 8 * (std::pow(D, 3)) * Rx0 * h * (
                             std::pow(nu, 3)) + 4 * (std::pow(D, 2)) * R0y * Rx0 * (std::pow(h, 2)) * nu) / (D * (2 * D + Rx0 * h) * (
            4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (std::pow(h, 2)))) - (
                     4 * D * Rx0 * h * nu) / (
                     4 * (std::pow(D, 2)) - 4 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 2 * D * R0y * h + 2 * D * Rx0 * h + R0y * Rx0 * (
                     std::pow(h, 2))) + 2);

    return {u10_c, u20_c, u30_c, u00_c, u11_c, u12_c, u21_c, u01_c};
}

std::vector<double> D11_coeffs_x(double R0y, double Rx0, double h, double D, double nu){
    return {(20 - (2 * D - Rx0 * h) / (2 * D + Rx0 * h) - (2 * D - R0y * h) / (2 * D + R0y * h)), - 8, 1,
            ((4 * D + 4 * D * nu) / (2 * D + Rx0 * h) - 8), ((4 * D + 4 * D * nu) / (2 * D + R0y * h) - 8), - 8,
            1, 2, (2 - (2 * D * nu) / (2 * D + Rx0 * h)),
            (2 - (2 * D * nu) / (2 * D + Rx0 * h) - (2 * D * nu) / (2 * D + R0y * h)),
            (2 - (2 * D * nu) / (2 * D + R0y * h))};
}

std::vector<double> D12_coeffs_x(double R0y, double h, double D, double nu){
    return {(20 - (2 * D - R0y * h) / (2 * D + R0y * h)), - 8, 1, - 8, 1,
            ((4 * D + 4 * D * nu) / (2 * D + R0y * h) - 8), - 8, 1, 2, 2, (2 - (2 * D * nu) / (2 * D + R0y * h)),
            (2 - (2 * D * nu) / (2 * D + R0y * h))};
}

std::vector<double> D20_coeffs_x(double Kx0, double Rx0, double h, double D, double nu){
    /*param Kx0:
    :param Rx0:
    :param h:
    :param D:
    :param nu:
    :return*/
    double u20_c = ((2 * Kx0 * Rx0 * (std::pow(h, 4)) + 4 * Kx0 * D * (std::pow(h, 3)) - 12 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 24 * (
            std::pow(D, 2)) * nu + 24 * (std::pow(D, 2))) / (D * (2 * D + Rx0 * h)) - (8 * D * nu) / (2 * D + Rx0 * h) - (
                     8 * (4 * D + 4 * D * nu)) / (2 * D + Rx0 * h) + 20);

    double u21_c = ((8 * (2 * D - Rx0 * h)) / (2 * D + Rx0 * h) + (8 * (std::pow(D, 2)) * nu - 24 * (std::pow(D, 2))) / (
            D * (2 * D + Rx0 * h)) - 8);

    double u22_c = ((2 * (std::pow(D, 2)) + Rx0 * h * D) / (D * (2 * D + Rx0 * h)) + 1);

    double u10_c = ((2 * (4 * D + 4 * D * nu)) / (2 * D + Rx0 * h) - (
            - 8 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 16 * (std::pow(D, 2)) * nu + 8 * (std::pow(D, 2))) / (D * (2 * D + Rx0 * h)) + (
                     16 * D * nu) / (2 * D + Rx0 * h) - 8);
    double u30_c = u10_c;

    double u40_c = ((- 2 * (std::pow(D, 2)) * (std::pow(nu, 2)) + 4 * (std::pow(D, 2)) * nu) / (D * (2 * D + Rx0 * h)) - (4 * D * nu) / (
            2 * D + Rx0 * h) + 1);
    double u00_c = u40_c;

    double u31_c = (2 - (4 * (std::pow(D, 2)) * nu - 8 * (std::pow(D, 2))) / (D * (2 * D + Rx0 * h)) - (2 * (2 * D - Rx0 * h)) / (
            2 * D + Rx0 * h));
    double u11_c = u31_c;
    return {u20_c, u21_c, u22_c, u10_c, u30_c, u40_c, u00_c, u31_c, u11_c};
}

std::vector<double> D21_coeffs_x(double Rx0, double h, double D, double nu){
    return {(20 - (2 * D - Rx0 * h) / (2 * D + Rx0 * h)), - 8, 1, + ((4 * D + 4 * D * nu) / (2 * D + Rx0 * h) - 8),
            - 8, - 8, 1, 1, 2, (2 - (2 * D * nu) / (2 * D + Rx0 * h)), (2 - (2 * D * nu) / (2 * D + Rx0 * h)), 2};
}

std::vector<double> D22_coeffs_x(){
    return {1, 2, -8, 2, 1, -8, 20, -8, 1, 2, -8, 2, 1};
}
#endif /* Dxx_coeffs_h */
