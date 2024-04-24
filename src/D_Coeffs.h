//
//  D_Coeffs.h
//  using_eigen
//
//  Created by Matthew Hamilton on 23/04/2024.
//

#ifndef D_Coeffs_h
#define D_Coeffs_h

#include <cmath>
#include <array>
#include <limits>

std::array<double,6>   D00_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,8>   D01_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,9>   D02_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,8>   D10_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,11>  D11_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,12>  D12_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,9>   D20_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,12>  D21_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);
std::array<double,13>  D22_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i);


std::array<double,6> D00_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    double u00_c = std::numeric_limits<double>::infinity();
    double u10_c = std::numeric_limits<double>::infinity();
    double u20_c = std::numeric_limits<double>::infinity();
    double u01_c = std::numeric_limits<double>::infinity();
    double u02_c = std::numeric_limits<double>::infinity();
    double u11_c = std::numeric_limits<double>::infinity();
    
    switch (i)
    {
        case 0:
            u00_c = ((2*K0y*Rx0*std::pow(R0y, 2)*std::pow(h, 6) + 4*K0y*D*std::pow(R0y, 2)*std::pow(h, 5) + 8*K0y*Rx0*D*R0y*std::pow(h, 5) - 8*K0y*std::pow(D, 2)*R0y*std::pow(h, 4)*std::pow(nu, 2) + 16*K0y*std::pow(D, 2)*R0y*std::pow(h, 4) - 12*Rx0*std::pow(D, 2)*R0y*std::pow(h, 2)*std::pow(nu, 2) + 24*Rx0*std::pow(D, 2)*R0y*std::pow(h, 2)*nu + 24*Rx0*std::pow(D, 2)*R0y*std::pow(h, 2) + 16*std::pow(D, 3)*R0y*h*std::pow(nu, 3) - 56*std::pow(D, 3)*R0y*h*std::pow(nu, 2) + 48*std::pow(D, 3)*R0y*h + 8*K0y*Rx0*std::pow(D, 2)*std::pow(h, 4) - 16*K0y*std::pow(D, 3)*std::pow(h, 3)*std::pow(nu, 2) + 16*K0y*std::pow(D, 3)*std::pow(h, 3) - 24*Rx0*std::pow(D, 3)*h*std::pow(nu, 2) + 48*Rx0*std::pow(D, 3)*h*nu + 48*Rx0*std::pow(D, 3)*h + 16*std::pow(D, 4)*std::pow(nu, 4) - 112*std::pow(D, 4)*std::pow(nu, 2) + 96*std::pow(D, 4))/(D*(2*D + R0y*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (8*(- 8*std::pow(D, 2)*std::pow(nu, 2) + 4*Rx0*h*D*nu + 8*std::pow(D, 2) + 4*Rx0*h*D))/(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)) - (4*D*nu)/(2*D + R0y*h) - (4*D*nu)/(2*D + Rx0*h) - (2*(8*std::pow(D, 2)*nu + 2*D*R0y*h*nu + 2*D*Rx0*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) - (8*(- 8*std::pow(D, 2)*std::pow(nu, 2) + 4*R0y*h*D*nu + 8*std::pow(D, 2) + 4*R0y*h*D))/(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)) + (2*Kx0*R0y*std::pow(Rx0, 2)*std::pow(h, 6) + 4*Kx0*D*std::pow(Rx0, 2)*std::pow(h, 5) + 8*Kx0*R0y*D*Rx0*std::pow(h, 5) - 8*Kx0*std::pow(D, 2)*Rx0*std::pow(h, 4)*std::pow(nu, 2) + 16*Kx0*std::pow(D, 2)*Rx0*std::pow(h, 4) - 12*R0y*std::pow(D, 2)*Rx0*std::pow(h, 2)*std::pow(nu, 2) + 24*R0y*std::pow(D, 2)*Rx0*std::pow(h, 2)*nu + 24*R0y*std::pow(D, 2)*Rx0*std::pow(h, 2) + 16*std::pow(D, 3)*Rx0*h*std::pow(nu, 3) - 56*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) + 48*std::pow(D, 3)*Rx0*h + 8*Kx0*R0y*std::pow(D, 2)*std::pow(h, 4) - 16*Kx0*std::pow(D, 3)*std::pow(h, 3)*std::pow(nu, 2) + 16*Kx0*std::pow(D, 3)*std::pow(h, 3) - 24*R0y*std::pow(D, 3)*h*std::pow(nu, 2) + 48*R0y*std::pow(D, 3)*h*nu + 48*R0y*std::pow(D, 3)*h + 16*std::pow(D, 4)*std::pow(nu, 4) - 112*std::pow(D, 4)*std::pow(nu, 2) + 96*std::pow(D, 4))/(D*(2*D + Rx0*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) + 20) ;
            break;
        case 1:
            u10_c =  ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) - (8*(4*std::pow(D, 2)*std::pow(nu, 2) - 4*std::pow(D, 2) + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)))/(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)) + (2*(8*std::pow(D, 2)*nu + 8*std::pow(D, 2) + 4*D*R0y*h + 4*D*R0y*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) + (32*std::pow(D, 4)*nu - 96*std::pow(D, 4) + 96*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) - 48*std::pow(D, 3)*R0y*h - 48*std::pow(D, 3)*Rx0*h + 16*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu - 24*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2) + 48*std::pow(D, 3)*R0y*h*std::pow(nu, 2) - 16*std::pow(D, 3)*R0y*h*std::pow(nu, 3) + 8*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (32*std::pow(D, 4)*nu + 64*std::pow(D, 4) - 96*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) + 32*std::pow(D, 4)*std::pow(nu, 4) + 32*std::pow(D, 3)*R0y*h + 32*std::pow(D, 3)*Rx0*h + 64*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu + 16*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2) - 32*std::pow(D, 3)*R0y*h*std::pow(nu, 2) - 16*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) - 16*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*std::pow(nu, 2) + 32*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) + (32*D*R0y*h*nu)/(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)) - 8) ;
            break;
        case 2:
            u20_c =  ((Rx0*D*std::pow(R0y, 2)*std::pow(h, 3) + 2*std::pow(D, 2)*std::pow(R0y, 2)*std::pow(h, 2) + 4*Rx0*std::pow(D, 2)*R0y*std::pow(h, 2) - 4*std::pow(D, 3)*R0y*h*std::pow(nu, 2) + 8*std::pow(D, 3)*R0y*h + 4*Rx0*std::pow(D, 3)*h - 8*std::pow(D, 4)*std::pow(nu, 2) + 8*std::pow(D, 4))/(D*(2*D + R0y*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (2*(4*nu*std::pow(D, 2) + 2*R0y*h*nu*D))/((2*D + R0y*h)*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + (32*std::pow(D, 4)*nu - 16*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) + 16*std::pow(D, 4)*std::pow(nu, 4) + 16*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu - 8*std::pow(D, 3)*R0y*h*std::pow(nu, 2) - 8*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) - 4*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*std::pow(nu, 2) + 8*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) + 1) ;
            break;
        case 3:
            u01_c = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (8*(4*std::pow(D, 2)*std::pow(nu, 2) - 4*std::pow(D, 2) - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)))/(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)) + (2*(8*std::pow(D, 2)*nu + 8*std::pow(D, 2) + 4*D*Rx0*h + 4*D*Rx0*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) + (32*std::pow(D, 4)*nu - 96*std::pow(D, 4) + 96*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) - 48*std::pow(D, 3)*R0y*h - 48*std::pow(D, 3)*Rx0*h + 16*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu - 24*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2) + 48*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) - 16*std::pow(D, 3)*Rx0*h*std::pow(nu, 3) + 8*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (32*std::pow(D, 4)*nu + 64*std::pow(D, 4) - 96*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) + 32*std::pow(D, 4)*std::pow(nu, 4) + 32*std::pow(D, 3)*R0y*h + 32*std::pow(D, 3)*Rx0*h + 16*std::pow(D, 3)*R0y*h*nu + 64*std::pow(D, 3)*Rx0*h*nu + 16*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2) - 16*std::pow(D, 3)*R0y*h*std::pow(nu, 2) - 32*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) - 16*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*std::pow(nu, 2) + 32*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) + (32*D*Rx0*h*nu)/(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2)) - 8) ;
            break;
        case 4:
            u02_c = ((R0y*D*std::pow(Rx0, 2)*std::pow(h, 3) + 2*std::pow(D, 2)*std::pow(Rx0, 2)*std::pow(h, 2) + 4*R0y*std::pow(D, 2)*Rx0*std::pow(h, 2) - 4*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) + 8*std::pow(D, 3)*Rx0*h + 4*R0y*std::pow(D, 3)*h - 8*std::pow(D, 4)*std::pow(nu, 2) + 8*std::pow(D, 4))/(D*(2*D + Rx0*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (2*(4*nu*std::pow(D, 2) + 2*Rx0*h*nu*D))/((2*D + R0y*h)*(2*D + Rx0*h)) - (4*D*nu)/(2*D + R0y*h) + (32*std::pow(D, 4)*nu - 16*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) + 16*std::pow(D, 4)*std::pow(nu, 4) + 16*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu - 8*std::pow(D, 3)*R0y*h*std::pow(nu, 2) - 8*std::pow(D, 3)*Rx0*h*std::pow(nu, 2) - 4*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*std::pow(nu, 2) + 8*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) + 1) ;
            break;
        case 5:
            u11_c = (2 - (2*(2*D - Rx0*h))/(2*D + Rx0*h) - (2*(12*std::pow(D, 2) + 2*D*R0y*h + 2*D*Rx0*h - R0y*Rx0*std::pow(h, 2)))/((2*D + R0y*h)*(2*D + Rx0*h)) - (32*std::pow(D, 4)*nu - 64*std::pow(D, 4) + 64*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) - 32*std::pow(D, 3)*R0y*h - 32*std::pow(D, 3)*Rx0*h + 16*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu - 16*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2) + 8*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (32*std::pow(D, 4)*nu - 64*std::pow(D, 4) + 64*std::pow(D, 4)*std::pow(nu, 2) - 32*std::pow(D, 4)*std::pow(nu, 3) - 32*std::pow(D, 3)*R0y*h - 32*std::pow(D, 3)*Rx0*h + 16*std::pow(D, 3)*R0y*h*nu + 16*std::pow(D, 3)*Rx0*h*nu - 16*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2) + 8*std::pow(D, 2)*R0y*Rx0*std::pow(h, 2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D, 2) - 4*std::pow(D, 2)*std::pow(nu, 2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h, 2))) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;
            break;
    }
    return {
        u00_c,u10_c,u20_c,u01_c,u02_c,u11_c
    };
}


std::array<double, 8> D01_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    double u01_c = std::numeric_limits<double>::infinity();
    double u11_c = std::numeric_limits<double>::infinity();
    double u21_c = std::numeric_limits<double>::infinity();
    double u00_c = std::numeric_limits<double>::infinity();
    double u02_c = std::numeric_limits<double>::infinity();
    double u03_c = std::numeric_limits<double>::infinity();
    double u12_c = std::numeric_limits<double>::infinity();
    double u10_c = std::numeric_limits<double>::infinity();
    
    switch (i) {
        case 0:
            u01_c =  ((4*std::pow(D,2)*std::pow(nu,2) - 4*std::pow(D,2) - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) - (4*D*nu)/(2*D + R0y*h) + (2*K0y*Rx0*std::pow(R0y,2)*std::pow(h,6) + 4*K0y*D*std::pow(R0y,2)*std::pow(h,5) + 8*K0y*Rx0*D*R0y*std::pow(h,5) - 8*K0y*std::pow(D,2)*R0y*std::pow(h,4)*std::pow(nu,2) + 16*K0y*std::pow(D,2)*R0y*std::pow(h,4) - 14*Rx0*std::pow(D,2)*R0y*std::pow(h,2)*std::pow(nu,2) + 28*Rx0*std::pow(D,2)*R0y*std::pow(h,2)*nu + 24*Rx0*std::pow(D,2)*R0y*std::pow(h,2) - 20*std::pow(D,3)*R0y*h*std::pow(nu,2) + 40*std::pow(D,3)*R0y*h*nu + 48*std::pow(D,3)*R0y*h + 8*K0y*Rx0*std::pow(D,2)*std::pow(h,4) - 16*K0y*std::pow(D,3)*std::pow(h,3)*std::pow(nu,2) + 16*K0y*std::pow(D,3)*std::pow(h,3) - 28*Rx0*std::pow(D,3)*h*std::pow(nu,2) + 56*Rx0*std::pow(D,3)*h*nu + 48*Rx0*std::pow(D,3)*h + 40*std::pow(D,4)*std::pow(nu,4) - 80*std::pow(D,4)*std::pow(nu,3) - 136*std::pow(D,4)*std::pow(nu,2) + 80*std::pow(D,4)*nu + 96*std::pow(D,4))/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (8*D*Rx0*h*nu)/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + 20) ;
            break;
        case 1:
            u11_c =  ((8*(2*D - R0y*h))/(2*D + R0y*h) + (32*std::pow(D,4)*nu - 96*std::pow(D,4) + 96*std::pow(D,4)*std::pow(nu,2) - 32*std::pow(D,4)*std::pow(nu,3) - 48*std::pow(D,3)*R0y*h - 48*std::pow(D,3)*Rx0*h + 16*std::pow(D,3)*R0y*h*nu + 16*std::pow(D,3)*Rx0*h*nu - 24*std::pow(D,2)*R0y*Rx0*std::pow(h,2) + 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - 8) ;
            break;
        case 2:
            u21_c =  ((Rx0*D*std::pow(R0y,2)*std::pow(h,3) + 2*std::pow(D,2)*std::pow(R0y,2)*std::pow(h,2) + 4*Rx0*std::pow(D,2)*R0y*std::pow(h,2) - 4*std::pow(D,3)*R0y*h*std::pow(nu,2) + 8*std::pow(D,3)*R0y*h + 4*Rx0*std::pow(D,3)*h - 8*std::pow(D,4)*std::pow(nu,2) + 8*std::pow(D,4))/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) + 1) ;
            break;
        case 3:
            u00_c =  ((- 8*std::pow(D,2)*std::pow(nu,2) + 4*R0y*h*D*nu + 8*std::pow(D,2) + 4*R0y*h*D)/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + (2*(- 8*std::pow(D,2)*std::pow(nu,2) + 4*Rx0*h*D*nu + 8*std::pow(D,2) + 4*Rx0*h*D))/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + (16*D*nu)/(2*D + R0y*h) - (32*std::pow(D,4)*nu + 32*std::pow(D,4) - 48*std::pow(D,4)*std::pow(nu,2) - 32*std::pow(D,4)*std::pow(nu,3) + 16*std::pow(D,4)*std::pow(nu,4) + 16*std::pow(D,3)*R0y*h + 16*std::pow(D,3)*Rx0*h + 16*std::pow(D,3)*R0y*h*nu + 32*std::pow(D,3)*Rx0*h*nu + 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) - 24*std::pow(D,3)*R0y*h*std::pow(nu,2) + 8*std::pow(D,3)*R0y*h*std::pow(nu,3) - 16*std::pow(D,3)*Rx0*h*std::pow(nu,2) - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*std::pow(nu,2) + 16*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - 8) ;
            break;
        case 4:
            u02_c =  ((2*(4*D + 4*D*nu))/(2*D + R0y*h) + (16*D*nu)/(2*D + R0y*h) - (64*std::pow(D,4)*nu + 32*std::pow(D,4) - 64*std::pow(D,4)*std::pow(nu,2) - 64*std::pow(D,4)*std::pow(nu,3) + 32*std::pow(D,4)*std::pow(nu,4) + 16*std::pow(D,3)*R0y*h + 16*std::pow(D,3)*Rx0*h + 32*std::pow(D,3)*R0y*h*nu + 32*std::pow(D,3)*Rx0*h*nu + 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) - 16*std::pow(D,3)*R0y*h*std::pow(nu,2) - 16*std::pow(D,3)*Rx0*h*std::pow(nu,2) - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*std::pow(nu,2) + 16*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - 8) ;
            break;
        case 5:
            u03_c =  ((16*std::pow(D,4)*nu - 8*std::pow(D,4)*std::pow(nu,2) - 16*std::pow(D,4)*std::pow(nu,3) + 8*std::pow(D,4)*std::pow(nu,4) + 8*std::pow(D,3)*R0y*h*nu + 8*std::pow(D,3)*Rx0*h*nu - 4*std::pow(D,3)*R0y*h*std::pow(nu,2) - 4*std::pow(D,3)*Rx0*h*std::pow(nu,2) - 2*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*std::pow(nu,2) + 4*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (4*D*nu)/(2*D + R0y*h) + 1) ;
            break;
        case 6:
            u12_c =  (2 - (16*std::pow(D,4)*nu - 32*std::pow(D,4) + 32*std::pow(D,4)*std::pow(nu,2) - 16*std::pow(D,4)*std::pow(nu,3) - 16*std::pow(D,3)*R0y*h - 16*std::pow(D,3)*Rx0*h + 8*std::pow(D,3)*R0y*h*nu + 8*std::pow(D,3)*Rx0*h*nu - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) + 4*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;
            break;
        case 7:
            u10_c =  ((2*(4*std::pow(D,2)*std::pow(nu,2) - 4*std::pow(D,2) + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)))/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) - (16*std::pow(D,4)*nu - 32*std::pow(D,4) + 32*std::pow(D,4)*std::pow(nu,2) - 16*std::pow(D,4)*std::pow(nu,3) - 16*std::pow(D,3)*R0y*h - 16*std::pow(D,3)*Rx0*h + 8*std::pow(D,3)*R0y*h*nu + 8*std::pow(D,3)*Rx0*h*nu - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) + 16*std::pow(D,3)*R0y*h*std::pow(nu,2) - 8*std::pow(D,3)*R0y*h*std::pow(nu,3) + 4*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + R0y*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (4*D*R0y*h*nu)/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + 2) ;
            break;
    }
    
    return {
        u01_c, u11_c, u21_c, u00_c, u02_c, u03_c, u12_c, u10_c
    };
}


std::array<double,9> D02_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    double u02_c = std::numeric_limits<double>::infinity();
    double u12_c = std::numeric_limits<double>::infinity();
    double u22_c = std::numeric_limits<double>::infinity();
    double u01_c = std::numeric_limits<double>::infinity();
    double u03_c = std::numeric_limits<double>::infinity();
    double u04_c = std::numeric_limits<double>::infinity();
    double u00_c = std::numeric_limits<double>::infinity();
    double u13_c = std::numeric_limits<double>::infinity();
    double u11_c = std::numeric_limits<double>::infinity();

    switch (i)
    {
        case 0:
            u02_c =  ((2*K0y*R0y*std::pow(h,4) + 4*K0y*D*std::pow(h,3) - 12*std::pow(D,2)*std::pow(nu,2) + 24*std::pow(D,2)*nu + 24*std::pow(D,2))/(D*(2*D + R0y*h)) - (8*D*nu)/(2*D + R0y*h) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) + 20) ;
            break;
        case 1:
            u12_c =  ((8*(2*D - R0y*h))/(2*D + R0y*h) + (8*std::pow(D,2)*nu - 24*std::pow(D,2))/(D*(2*D + R0y*h)) - 8) ;
            break;
        case 2:
            u22_c =  ((2*std::pow(D,2) + R0y*h*D)/(D*(2*D + R0y*h)) + 1) ;
            break;
        case 3:
            u01_c =  ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (- 8*std::pow(D,2)*std::pow(nu,2) + 16*std::pow(D,2)*nu + 8*std::pow(D,2))/(D*(2*D + R0y*h)) + (16*D*nu)/(2*D + R0y*h) - 8) ;
            break;
        case 4:
            u03_c =  ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (- 8*std::pow(D,2)*std::pow(nu,2) + 16*std::pow(D,2)*nu + 8*std::pow(D,2))/(D*(2*D + R0y*h)) + (16*D*nu)/(2*D + R0y*h) - 8) ;;
            break;
        case 5:
            u04_c =  ((- 2*std::pow(D,2)*std::pow(nu,2) + 4*std::pow(D,2)*nu)/(D*(2*D + R0y*h)) - (4*D*nu)/(2*D + R0y*h) + 1) ;
            break;
        case 6:
            u00_c =  ((- 2*std::pow(D,2)*std::pow(nu,2) + 4*std::pow(D,2)*nu)/(D*(2*D + R0y*h)) - (4*D*nu)/(2*D + R0y*h) + 1);
            break;
        case 7:
            u13_c =  (2 - (4*std::pow(D,2)*nu - 8*std::pow(D,2))/(D*(2*D + R0y*h)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;
            break;
        case 8:
            u11_c =  (2 - (4*std::pow(D,2)*nu - 8*std::pow(D,2))/(D*(2*D + R0y*h)) - (2*(2*D - R0y*h))/(2*D + R0y*h));
            break;
    }

    return {
        u02_c, u12_c, u22_c, u01_c, u03_c, u04_c, u00_c, u13_c, u11_c
    };
}

std::array<double,8> D10_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    double u10_c = std::numeric_limits<double>::infinity();
    double u20_c = std::numeric_limits<double>::infinity();
    double u30_c = std::numeric_limits<double>::infinity();
    double u00_c = std::numeric_limits<double>::infinity();
    double u11_c = std::numeric_limits<double>::infinity();
    double u12_c = std::numeric_limits<double>::infinity();
    double u21_c = std::numeric_limits<double>::infinity();
    double u01_c = std::numeric_limits<double>::infinity();

    switch (i)
    {
        case 0:
            u10_c = ((4*std::pow(D,2)*std::pow(nu,2) - 4*std::pow(D,2) + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) - (4*D*nu)/(2*D + Rx0*h) + (2*Kx0*R0y*std::pow(Rx0,2)*std::pow(h,6) + 4*Kx0*D*std::pow(Rx0,2)*std::pow(h,5) + 8*Kx0*R0y*D*Rx0*std::pow(h,5) - 8*Kx0*std::pow(D,2)*Rx0*std::pow(h,4)*std::pow(nu,2) + 16*Kx0*std::pow(D,2)*Rx0*std::pow(h,4) - 14*R0y*std::pow(D,2)*Rx0*std::pow(h,2)*std::pow(nu,2) + 28*R0y*std::pow(D,2)*Rx0*std::pow(h,2)*nu + 24*R0y*std::pow(D,2)*Rx0*std::pow(h,2) - 20*std::pow(D,3)*Rx0*h*std::pow(nu,2) + 40*std::pow(D,3)*Rx0*h*nu + 48*std::pow(D,3)*Rx0*h + 8*Kx0*R0y*std::pow(D,2)*std::pow(h,4) - 16*Kx0*std::pow(D,3)*std::pow(h,3)*std::pow(nu,2) + 16*Kx0*std::pow(D,3)*std::pow(h,3) - 28*R0y*std::pow(D,3)*h*std::pow(nu,2) + 56*R0y*std::pow(D,3)*h*nu + 48*R0y*std::pow(D,3)*h + 40*std::pow(D,4)*std::pow(nu,4) - 80*std::pow(D,4)*std::pow(nu,3) - 136*std::pow(D,4)*std::pow(nu,2) + 80*std::pow(D,4)*nu + 96*std::pow(D,4))/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (8*D*R0y*h*nu)/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + 20) ;
            break;
        case 1:
            u20_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) + (16*D*nu)/(2*D + Rx0*h) - (64*std::pow(D,4)*nu + 32*std::pow(D,4) - 64*std::pow(D,4)*std::pow(nu,2) - 64*std::pow(D,4)*std::pow(nu,3) + 32*std::pow(D,4)*std::pow(nu,4) + 16*std::pow(D,3)*R0y*h + 16*std::pow(D,3)*Rx0*h + 32*std::pow(D,3)*R0y*h*nu + 32*std::pow(D,3)*Rx0*h*nu + 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) - 16*std::pow(D,3)*R0y*h*std::pow(nu,2) - 16*std::pow(D,3)*Rx0*h*std::pow(nu,2) - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*std::pow(nu,2) + 16*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - 8) ;
            break;
        case 2:
            u30_c = ((16*std::pow(D,4)*nu - 8*std::pow(D,4)*std::pow(nu,2) - 16*std::pow(D,4)*std::pow(nu,3) + 8*std::pow(D,4)*std::pow(nu,4) + 8*std::pow(D,3)*R0y*h*nu + 8*std::pow(D,3)*Rx0*h*nu - 4*std::pow(D,3)*R0y*h*std::pow(nu,2) - 4*std::pow(D,3)*Rx0*h*std::pow(nu,2) - 2*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*std::pow(nu,2) + 4*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (4*D*nu)/(2*D + Rx0*h) + 1) ;
            break;
        case 3:
            u00_c = ((2*(- 8*std::pow(D,2)*std::pow(nu,2) + 4*R0y*h*D*nu + 8*std::pow(D,2) + 4*R0y*h*D))/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + (- 8*std::pow(D,2)*std::pow(nu,2) + 4*Rx0*h*D*nu + 8*std::pow(D,2) + 4*Rx0*h*D)/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + (16*D*nu)/(2*D + Rx0*h) - (32*std::pow(D,4)*nu + 32*std::pow(D,4) - 48*std::pow(D,4)*std::pow(nu,2) - 32*std::pow(D,4)*std::pow(nu,3) + 16*std::pow(D,4)*std::pow(nu,4) + 16*std::pow(D,3)*R0y*h + 16*std::pow(D,3)*Rx0*h + 32*std::pow(D,3)*R0y*h*nu + 16*std::pow(D,3)*Rx0*h*nu + 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) - 16*std::pow(D,3)*R0y*h*std::pow(nu,2) - 24*std::pow(D,3)*Rx0*h*std::pow(nu,2) + 8*std::pow(D,3)*Rx0*h*std::pow(nu,3) - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*std::pow(nu,2) + 16*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - 8) ;
            break;
        case 4:
            u11_c = ((8*(2*D - Rx0*h))/(2*D + Rx0*h) + (32*std::pow(D,4)*nu - 96*std::pow(D,4) + 96*std::pow(D,4)*std::pow(nu,2) - 32*std::pow(D,4)*std::pow(nu,3) - 48*std::pow(D,3)*R0y*h - 48*std::pow(D,3)*Rx0*h + 16*std::pow(D,3)*R0y*h*nu + 16*std::pow(D,3)*Rx0*h*nu - 24*std::pow(D,2)*R0y*Rx0*std::pow(h,2) + 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - 8) ;
            break;
        case 5:
            u12_c = ((R0y*D*std::pow(Rx0,2)*std::pow(h,3) + 2*std::pow(D,2)*std::pow(Rx0,2)*std::pow(h,2) + 4*R0y*std::pow(D,2)*Rx0*std::pow(h,2) - 4*std::pow(D,3)*Rx0*h*std::pow(nu,2) + 8*std::pow(D,3)*Rx0*h + 4*R0y*std::pow(D,3)*h - 8*std::pow(D,4)*std::pow(nu,2) + 8*std::pow(D,4))/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) + 1) ;
            break;
        case 6:
            u21_c = (2 - (16*std::pow(D,4)*nu - 32*std::pow(D,4) + 32*std::pow(D,4)*std::pow(nu,2) - 16*std::pow(D,4)*std::pow(nu,3) - 16*std::pow(D,3)*R0y*h - 16*std::pow(D,3)*Rx0*h + 8*std::pow(D,3)*R0y*h*nu + 8*std::pow(D,3)*Rx0*h*nu - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) + 4*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (2*(2*D - Rx0*h))/(2*D + Rx0*h)) ;
            break;
        case 7:
            u01_c = ((2*(4*std::pow(D,2)*std::pow(nu,2) - 4*std::pow(D,2) - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)))/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) - (16*std::pow(D,4)*nu - 32*std::pow(D,4) + 32*std::pow(D,4)*std::pow(nu,2) - 16*std::pow(D,4)*std::pow(nu,3) - 16*std::pow(D,3)*R0y*h - 16*std::pow(D,3)*Rx0*h + 8*std::pow(D,3)*R0y*h*nu + 8*std::pow(D,3)*Rx0*h*nu - 8*std::pow(D,2)*R0y*Rx0*std::pow(h,2) + 16*std::pow(D,3)*Rx0*h*std::pow(nu,2) - 8*std::pow(D,3)*Rx0*h*std::pow(nu,3) + 4*std::pow(D,2)*R0y*Rx0*std::pow(h,2)*nu)/(D*(2*D + Rx0*h)*(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2))) - (4*D*Rx0*h*nu)/(4*std::pow(D,2) - 4*std::pow(D,2)*std::pow(nu,2) + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*std::pow(h,2)) + 2) ;
            break;
    }
    
    return {
        u10_c, u20_c, u30_c, u00_c, u11_c, u12_c, u21_c, u01_c
    };
}

std::array<double,11> D11_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    return {
        (20 - (2*D - Rx0*h)/(2*D + Rx0*h) - (2*D - R0y*h)/(2*D + R0y*h)), 
        - 8,
        1,
        ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8), 
        ((4*D + 4*D*nu)/(2*D + R0y*h) - 8),
        - 8,
        1, 
        2,
        (2 - (2*D*nu)/(2*D + Rx0*h)),
        (2 - (2*D*nu)/(2*D + Rx0*h) - (2*D*nu)/(2*D + R0y*h)),
        (2 - (2*D*nu)/(2*D + R0y*h))
    };    
}

std::array<double,12> D12_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    return {
        (20 - (2*D - R0y*h)/(2*D + R0y*h)),
        - 8,
        1,
        - 8,
        1,
        ((4*D + 4*D*nu)/(2*D + R0y*h) - 8),
        - 8,
        1,
        2,
        2,
        (2 - (2*D*nu)/(2*D + R0y*h)),
        (2 - (2*D*nu)/(2*D + R0y*h))
    };
}

std::array<double,9>  D20_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    double u20_c = std::numeric_limits<double>::infinity();
    double u21_c = std::numeric_limits<double>::infinity();
    double u22_c = std::numeric_limits<double>::infinity();
    double u10_c = std::numeric_limits<double>::infinity();
    double u30_c = std::numeric_limits<double>::infinity();
    double u40_c = std::numeric_limits<double>::infinity();
    double u00_c = std::numeric_limits<double>::infinity();
    double u31_c = std::numeric_limits<double>::infinity();
    double u11_c = std::numeric_limits<double>::infinity();
    
    switch (i)
    {
        case 0:
            u20_c = ((2*Kx0*Rx0*std::pow(h,4) + 4*Kx0*D*std::pow(h,3) - 12*std::pow(D,2)*std::pow(nu,2) + 24*std::pow(D,2)*nu + 24*std::pow(D,2))/(D*(2*D + Rx0*h)) - (8*D*nu)/(2*D + Rx0*h) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) + 20) ;
            break;
        case 1:
            u21_c = ((8*(2*D - Rx0*h))/(2*D + Rx0*h) + (8*std::pow(D,2)*nu - 24*std::pow(D,2))/(D*(2*D + Rx0*h)) - 8) ;
            break;
        case 2:
            u22_c = ((2*std::pow(D,2) + Rx0*h*D)/(D*(2*D + Rx0*h)) + 1) ;
            break;
        case 3:
            u10_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) - (- 8*std::pow(D,2)*std::pow(nu,2) + 16*std::pow(D,2)*nu + 8*std::pow(D,2))/(D*(2*D + Rx0*h)) + (16*D*nu)/(2*D + Rx0*h) - 8) ;
            break;
        case 4:
            u30_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) - (- 8*std::pow(D,2)*std::pow(nu,2) + 16*std::pow(D,2)*nu + 8*std::pow(D,2))/(D*(2*D + Rx0*h)) + (16*D*nu)/(2*D + Rx0*h) - 8) ;;
            break;
        case 5:
            u40_c = ((- 2*std::pow(D,2)*std::pow(nu,2) + 4*std::pow(D,2)*nu)/(D*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + 1) ;
            break;
        case 6:
            u00_c = ((- 2*std::pow(D,2)*std::pow(nu,2) + 4*std::pow(D,2)*nu)/(D*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + 1) ;;
            break;
        case 7:
            u31_c = (2 - (4*std::pow(D,2)*nu - 8*std::pow(D,2))/(D*(2*D + Rx0*h)) - (2*(2*D - Rx0*h))/(2*D + Rx0*h)) ;
            break;
        case 8:
            u11_c = (2 - (4*std::pow(D,2)*nu - 8*std::pow(D,2))/(D*(2*D + Rx0*h)) - (2*(2*D - Rx0*h))/(2*D + Rx0*h));
            break;
    }
    
    return {
        u20_c, u21_c, u22_c, u10_c, u30_c, u40_c, u00_c, u31_c, u11_c
    };
}

std::array<double,12> D21_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    return {
        (20 - (2*D - Rx0*h)/(2*D + Rx0*h)),
        -8,
        1,
        + ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8),
        - 8,
        - 8,
        1,
        1,
        2,
        (2 - (2*D*nu)/(2*D + Rx0*h)),
        (2 - (2*D*nu)/(2*D + Rx0*h)),
        2
    };
}

std::array<double, 13> D22_coeffs(double K0y, double R0y, double Kx0, double Rx0, double h, double D, double nu, int i)
{
    return {1, 2, -8, 2, 1, -8, 20, -8, 1, 2, -8, 2, 1};
}

#endif /* D_Coeffs_h */
