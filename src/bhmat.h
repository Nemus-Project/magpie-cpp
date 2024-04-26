//
//  bhmat.h
//
//
//  Created by Matthew Hamilton on 24/04/2024.
//

#ifndef bhmat_h
#define bhmat_h

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

#include "range.h"
#include "spdiags.h"
#include "D_Coeffs.h"


void bhmat(Eigen::SparseMatrix<double>& biharm,
           std::array<double,8>BCs,
           std::array<int,2> Nxy,
           double h, double Lz, double E, double nu)
{
    const auto&  [K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = BCs;
    const int Nx = Nxy[0];
    const int Ny = Nxy[1];
    const double D = E * (Lz*Lz*Lz) / 12 / (1-(nu*nu));
    
    Eigen::VectorXd a0(Ny);
    Eigen::VectorXd a1(Ny-1);
    Eigen::VectorXd a2(Ny-2);
    
    //// dm2Ny  // pad zeros at the end
    double D20u00 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,6)[6];
    double D21u01 = D21_coeffs(Rx0,0,0,0,h,D,nu,7)[7];
    double D22u02 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[4];
    double D2Nu0N = D20_coeffs(KxL,RxL,0,0,h,D,nu,6)[6];
    double D2Nm1u0Nm1 = D21_coeffs(RxL,0,0,0,h,D,nu,7)[7];
    
    Eigen::VectorXd dm2Ny0 = D22u02*a0;
    dm2Ny0(0) = D20u00;
    dm2Ny0(1) = D21u01;
    dm2Ny0(Ny-2) = D2Nm1u0Nm1;
    dm2Ny0(Ny-1) = D2Nu0N;
    
    double D10u30 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,2)[2];
    double D11u31 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,6)[6];
    double D12u32 = D12_coeffs(RLy,0,0,0,h,D,nu,7)[7];
    double D1Nu3N = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,2)[2];
    double D1Nm1u3Nm1 = D11_coeffs(RLy,RxL,0,0,h,D,nu,6)[6];
    
    Eigen::VectorXd dm2Ny1 = D12u32*a0;
    dm2Ny1(0) = D10u30;
    dm2Ny1(1) = D11u31;
    dm2Ny1(Ny-2) = D1Nm1u3Nm1;
    dm2Ny1(Ny-1) = D1Nu3N;
    
    double D00u20 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,2)[2];
    double D01u21 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,2)[2];
    double D02u22 = D02_coeffs(KLy,RLy,0,0,h,D,nu,2)[2];
    double D0Nu2N = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,2)[2];
    double D0Nm1u2Nm1 = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,2)[2];
    
    Eigen::VectorXd dm2Ny2 = D02u22*a0;
    dm2Ny2(0) = D00u20;
    dm2Ny2(1) = D01u21;
    dm2Ny2(Ny-2) = D0Nm1u2Nm1;
    dm2Ny2(Ny-1) = D0Nu2N;
    
    Eigen::VectorXd Dm2Ny(Ny*Nx);
    Dm2Ny << dm2Ny0.replicate(Nx-3,1),dm2Ny1,dm2Ny2,Eigen::VectorXd(Ny);
    
    //// dmNym1 // pad zeros at the end
    double D11u00 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,9)[9];
    double D12u01 = D12_coeffs(R0y,0,0,0,h,D,nu,10)[10];
    double D1Nu0Nm1 = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,7)[7];
    double D1Nm1u0Nm2 = D11_coeffs(R0y,RxL,0,0,h,D,nu,10)[10];
    
    
    Eigen::VectorXd dmNym10 = D12u01*a1;
    dmNym10(0) = D11u00;
    dmNym10(Ny-1-2) = D1Nm1u0Nm2;
    dmNym10(Ny-2) = D1Nu0Nm1;
    
    double D21u10 = D21_coeffs(Rx0,0,0,0,h,D,nu,10)[10];
    double D22u11 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[1];
    double D2Nu1Nm1 = D20_coeffs(KxL,RxL,0,0,h,D,nu,8)[8];
    double D2Nm1u1Nm2 = D21_coeffs(RxL,0,0,0,h,D,nu,11)[11];
    
    Eigen::VectorXd dmNym11 = D22u11*a0;
    dmNym11(0) = D21u10;
    dmNym11(Ny-1-2) = D2Nm1u1Nm2;
    dmNym11(Ny-2) = D2Nu1Nm1;
    
    double D11u20 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,8)[8];
    double D12u21 = D12_coeffs(RLy,0,0,0,h,D,nu,9)[9];
    double D1Nm1u2Nm2 = D11_coeffs(RLy,RxL,0,0,h,D,nu,7)[7];
    double D1Nu2Nm1 = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,6)[6];
    
    Eigen::VectorXd dmNym1M1 = D12u21*a1;
    dmNym1M1(0) = D11u20;
    dmNym1M1(Ny-1-2) = D1Nm1u2Nm2;
    dmNym1M1(Ny-2) = D1Nu2Nm1;
    
    double D01u10 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,7)[7];
    double D02u11 = D02_coeffs(KLy,RLy,0,0,h,D,nu,8)[8];
    double D0Nu1Nm1 = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,5)[5];
    double D0Nm1u1Nm2 = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,6)[6];
    
    Eigen::VectorXd dmNym1M = D02u11*a1;
    dmNym1M(0) = D01u10;
    dmNym1M(Ny-1-2) = D0Nm1u1Nm2;
    dmNym1M(Ny-2) = D0Nu1Nm1;
    
    
    Eigen::VectorXd DmNym1(Ny*Nx);
    DmNym1 << dmNym10,0,dmNym11.replicate(Nx-3,1),dmNym1M1,0,dmNym1M,0;
    //// dmNy   // pad zeros at the end
    
    double D10u00 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,3)[3];
    double D11u01 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,4)[4];
    double D12u02 = D12_coeffs(R0y,0,0,0,h,D,nu,5)[5];
    double D1Nu0N = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,3)[3];
    double D1Nm1u0Nm1 = D11_coeffs(R0y,RxL,0,0,h,D,nu,4)[4];
    
    Eigen::VectorXd dmNy0 = D12u02*a0;
    dmNy0(0) = D10u00;
    dmNy0(1) = D11u01;
    dmNy0(Ny-2) = D1Nm1u0Nm1;
    dmNy0(Ny-1) = D1Nu0N;
    
    double D20u10 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,3)[3];
    double D21u11 = D21_coeffs(Rx0,0,0,0,h,D,nu,4)[4];
    double D22u12 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[5];
    double D2Nu1N = D20_coeffs(KxL,RxL,0,0,h,D,nu,3)[3];
    double D2Nm1u1Nm1 = D21_coeffs(RxL,0,0,0,h,D,nu,4)[4];
    
    Eigen::VectorXd dmNy1 = D22u12*a0;
    dmNy1(0) = D20u10;
    dmNy1(1) = D21u11;
    dmNy1(Ny-2) = D2Nm1u1Nm1;
    dmNy1(Ny-1) = D2Nu1N;
    
    double D10u20 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,1)[1];
    double D11u21 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,5)[5];
    double D12u22 = D12_coeffs(RLy,0,0,0,h,D,nu,6)[6];
    double D1Nu2N = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,1)[1];
    double D1Nm1u2Nm1 = D11_coeffs(RLy,RxL,0,0,h,D,nu,5)[5];
    
    Eigen::VectorXd dmNyM1 = D12u22*a0;
    dmNyM1(0) = D10u20;
    dmNyM1(1) = D11u21;
    dmNyM1(Ny-2) = D1Nm1u2Nm1;
    dmNyM1(Ny-1) = D1Nu2N;
    
    double D00u10 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,1)[1];
    double D01u11 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,1)[1];
    double D02u12 = D02_coeffs(KLy,RLy,0,0,h,D,nu,1)[1];
    double D0Nu1N = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,1)[1];
    double D0Nm1u1Nm1 = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,1)[1];
    
    Eigen::VectorXd dmNyM = D02u12*a0;
    dmNyM(0) = D00u10;
    dmNyM(1) = D01u11;
    dmNyM(Ny-2) = D0Nm1u1Nm1;
    dmNyM(Ny-1) = D0Nu1N;
    
    Eigen::VectorXd DmNy(Ny*Nx);
    DmNy << dmNy0,dmNy1.replicate((Nx-3),1),dmNyM1,dmNyM;
    ////assert(all(abs(diag(biHarm,Ny+1) - DmNy) <= eps), "D0 incorrect");
    
    
    //// dmNyp1 // pad zeros at the
    double D10u01 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,7)[7];
    double D11u02 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,10)[10];
    double D12u03 = D12_coeffs(R0y,0,0,0,h,D,nu,11)[11];
    double D1Nm1u0N = D11_coeffs(R0y,RxL,0,0,h,D,nu,9)[9];
    
    Eigen::VectorXd dmNy10 = D12u03*a1;
    dmNy10(0) = D10u01;
    dmNy10(1) = D11u02;
    dmNy10(Ny-2) = D1Nm1u0N;
    
    double D20u11 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,8)[8];
    double D21u12 = D21_coeffs(Rx0,0,0,0,h,D,nu,11)[11];
    double D22u13 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[9];
    double D2Nm1u1N = D21_coeffs(RxL,0,0,0,h,D,nu,10)[10];
    
    
    Eigen::VectorXd dmNy11 = D22u13*a0;
    dmNy11(1) = D20u11;
    dmNy11(2) = D21u12;
    dmNy11(Ny-1) = D2Nm1u1N;
    
    double D10u21 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,6)[6];
    double D11u22 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,7)[7];
    double D12u23 = D12_coeffs(RLy,0,0,0,h,D,nu,8)[8];
    double D1Nm1u2N = D11_coeffs(RLy,RxL,0,0,h,D,nu,8)[8];
    
    
    Eigen::VectorXd dmNy1M1 = D12u23*a1;
    dmNy1M1(0) = D10u21;
    dmNy1M1(1) = D11u22;
    dmNy1M1(Ny-2) = D1Nm1u2N;
    
    double D00u11 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,5)[5];
    double D01u12 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,6)[6];
    double D02u13 = D02_coeffs(KLy,RLy,0,0,h,D,nu,7)[7];
    double D0Nm1u1N = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,7)[7];
    
    Eigen::VectorXd dmNy1M = D02u13*a1;
    dmNy1M(0) = D00u11;
    dmNy1M(1) = D01u12;
    dmNy1M(Ny-2) = D0Nm1u1N;
    
    
    Eigen::VectorXd DmNy1(Ny*Nx);
    DmNy1 << 0,dmNy10,dmNy11.replicate(Nx-3,1),0,dmNy1M1,0,dmNy1M;
    //// dm2   // pad zeros at the end
    double D02u00 = D02_coeffs(K0y,R0y,0,0,h,D,nu,6)[6];
    double D0Nu0Nm2 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,4)[4];
    double D0Nm1u0Nm3 = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,5)[5];
    
    
    Eigen::VectorXd dm20 = D02u00*a2;
    dm20(Ny-2-2) = D0Nm1u0Nm3;
    dm20(Ny-1-2) = D0Nu0Nm2;
    
    double D12u10 = D12_coeffs(R0y,0,0,0,h,D,nu,4)[4];
    double D1Nu1Nm2 = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,5)[5];
    double D1Nm1u1Nm3 = D11_coeffs(R0y,RxL,0,0,h,D,nu,2)[2];
    
    Eigen::VectorXd dm21 = D12u10*a2 ;
    dm21(Ny-2-2) = D1Nm1u1Nm3;
    dm21(Ny-1-2) = D1Nu1Nm2;
    
    double D22u20 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[0];
    double D2Nu2Nm2 = D20_coeffs(KxL,RxL,0,0,h,D,nu,2)[2];
    double D2Nm1u2Nm3 = D21_coeffs(RxL,0,0,0,h,D,nu,2)[2];
    
    
    Eigen::VectorXd dm22 = D22u20*a0;
    dm22(Ny-2-2) = D2Nm1u2Nm3;
    dm22(Ny-1-2) = D2Nu2Nm2;
    
    D12u10 = D12_coeffs(RLy,0,0,0,h,D,nu,4)[4];
    D1Nu1Nm2 = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,5)[5];
    D1Nm1u1Nm3 = D11_coeffs(RLy,RxL,0,0,h,D,nu,2)[2];
    
    
    Eigen::VectorXd dm2M1 = D12u10*a2 ;
    dm2M1(Ny-2-2) = D1Nm1u1Nm3;
    dm2M1(Ny-1-2) = D1Nu1Nm2;
    
    D02u00 = D02_coeffs(KLy,RLy,0,0,h,D,nu,6)[6];
    D0Nu0Nm2 = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,4)[4];
    D0Nm1u0Nm3 = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,5)[5];
    
    Eigen::VectorXd dm2M = D02u00*a2;
    dm2M(Ny-2-2) = D0Nm1u0Nm3;
    dm2M(Ny-1-2) = D0Nu0Nm2;
    
    
    Eigen::VectorXd Dm2(Ny*Nx);
    Dm2 << dm20,0,0,dm21,0,0,dm22.replicate(Nx-4,1),dm2M1,0,0,dm2M,0,0;
    //// dm1   // pad zeros at the end
    
    
    double D01u00 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,3)[3];
    double D02u01 = D02_coeffs(K0y,R0y,0,0,h,D,nu,3)[3];
    double D0Nu0Nm1 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,3)[3];
    double D0Nm1u0Nm2 = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,4)[4];
    
    Eigen::VectorXd dm10 = D02u01*a1;
    dm10(0) = D01u00;
    dm10(Ny-1-2) = D0Nm1u0Nm2;
    dm10(Ny-2) = D0Nu0Nm1;
    
    double D11u10 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,3)[3];
    double D12u11 = D12_coeffs(R0y,0,0,0,h,D,nu,3)[3];
    double D1Nu1Nm1 = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,4)[4];
    double D1Nm1u1Nm2 = D11_coeffs(R0y,RxL,0,0,h,D,nu,1)[1];
    
    
    Eigen::VectorXd dm11 = D12u11*a1;
    dm11(0) = D11u10;
    dm11(Ny-1-2) = D1Nm1u1Nm2;
    dm11(Ny-2) = D1Nu1Nm1;
    
    double D21u20 = D21_coeffs(Rx0,0,0,0,h,D,nu,3)[3];
    double D22u21 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[2];
    double D2Nu2Nm1 = D20_coeffs(KxL,RxL,0,0,h,D,nu,1)[1];
    double D2Nm1u2Nm2 = D21_coeffs(RxL,0,0,0,h,D,nu,1)[1];
    
    
    Eigen::VectorXd dm12 = D22u21*a0;
    dm12(0) = D21u20;
    dm12(Ny-1-2) = D2Nm1u2Nm2;
    dm12(Ny-2) = D2Nu2Nm1;
    
    D11u10 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,3)[3];
    D12u11 = D12_coeffs(RLy,0,0,0,h,D,nu,3)[3];
    D1Nu1Nm1 = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,4)[4];
    D1Nm1u1Nm2 = D11_coeffs(RLy,RxL,0,0,h,D,nu,1)[1];
    
    
    Eigen::VectorXd dm1M1 = D12u11*a1;
    dm1M1(0) = D11u10;
    dm1M1(Ny-1-2) = D1Nm1u1Nm2;
    dm1M1(Ny-2) = D1Nu1Nm1;
    
    D01u00 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,3)[3];
    D02u01 = D02_coeffs(KLy,RLy,0,0,h,D,nu,3)[3];
    D0Nu0Nm1 = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,3)[3];
    D0Nm1u0Nm2 = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,4)[4];

    
    Eigen::VectorXd dm1M = D02u01*a1;
    dm1M(0) = D01u00;
    dm1M(Ny-1-2) = D0Nm1u0Nm2;
    dm1M(Ny-2) = D0Nu0Nm1;
    
    
    Eigen::VectorXd Dm1(Ny*Nx);
    Dm1 << dm10,0,dm11,0,dm12.replicate(Nx-4,1),dm1M1,0,dm1M,0;
    //// d00
    
    double D00u00 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,0)[0];
    double D01u01 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,0)[0];
    double D02u02 = D02_coeffs(K0y,R0y,0,0,h,D,nu,0)[0];
    double D0Nu0N = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,0)[0];
    double D0Nm1u0Nm1 = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,0)[0];
    
    Eigen::VectorXd d00 = D02u02*a0;
    d00(0) = D00u00;
    d00(1) = D01u01;
    d00(Ny-2) = D0Nm1u0Nm1;
    d00(Ny-1) = D0Nu0N;
    
    double D10u10 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,0)[0];
    double D11u11 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,0)[0];
    double D12u12 = D12_coeffs(R0y,0,0,0,h,D,nu,0)[0];
    double D1Nu1N = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,0)[0];
    double D1Nm1u1Nm1 = D11_coeffs(R0y,RxL,0,0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d01 = D12u12*a0;
    d01(0) = D10u10;
    d01(1) = D11u11;
    d01(Ny-2) = D1Nm1u1Nm1;
    d01(Ny-1) = D1Nu1N;
    
    double D20u20 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,0)[0];
    double D21u21 = D21_coeffs(Rx0,0,0,0,h,D,nu,0)[0];
    double D22u22 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[6];
    double D2Nu2N = D20_coeffs(KxL,RxL,0,0,h,D,nu,0)[0];
    double D2Nm1u2Nm1 = D21_coeffs(RxL,0,0,0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d02 = D22u22*a0;
    d02(0) = D20u20;
    d02(1) = D21u21;
    d02(Ny-2) = D2Nm1u2Nm1;
    d02(Ny-1) = D2Nu2N;
    
    
    D10u10 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,0)[0];
    D11u11 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,0)[0];
    D12u12 = D12_coeffs(RLy,0,0,0,h,D,nu,0)[0];
    D1Nu1N = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,0)[0];
    D1Nm1u1Nm1 = D11_coeffs(RLy,RxL,0,0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d0Mm = D12u12*a0;
    d0Mm(0) = D10u10;
    d0Mm(1) = D11u11;
    d0Mm(Ny-2) = D1Nm1u1Nm1;
    d0Mm(Ny-1) = D1Nu1N;
    
    D00u00 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,0)[0];
    D01u01 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,0)[0];
    D02u02 = D02_coeffs(KLy,RLy,0,0,h,D,nu,0)[0];
    D0Nu0N = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,0)[0];
    D0Nm1u0Nm1 = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d0M = D02u02*a0;
    d0M(0) = D00u00;
    d0M(1) = D01u01;
    d0M(Ny-2) = D0Nm1u0Nm1;
    d0M(Ny-1) = D0Nu0N;
    
    Eigen::VectorXd D0(Ny*Nx);
    D0 << d00,d01,d02.replicate((Nx-4),1),d0Mm,d0M;
    // //assert(all(abs(diag(biHarm,0) - D0) <= eps), "D0 incorrect");
    
    
    //// dp1   // pad zeros at the start
    double D00u01 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,3)[3];
    double D01u02 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,4)[4];
    double D02u03 = D02_coeffs(K0y,R0y,0,0,h,D,nu,4)[4];
    double D0Nm1u0N = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,3)[3];
    
    Eigen::VectorXd d10 = D02u03*a1;
    d10(0) = D00u01;
    d10(1) = D01u02;
    d10(Ny-2) = D0Nm1u0N;
    
    double D10u11 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,4)[4];
    double D11u12 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,1)[1];
    double D12u13 = D12_coeffs(R0y,0,0,0,h,D,nu,1)[1];
    double D1Nm1u1N = D11_coeffs(R0y,RxL,0,0,h,D,nu,3)[3];
    
    Eigen::VectorXd d11 = D12u13*a1;
    d11(0) = D10u11;
    d11(1) = D11u12;
    d11(Ny-2) = D1Nm1u1N;
    
    double D20u21 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,1)[1];
    double D21u22 = D21_coeffs(Rx0,0,0,0,h,D,nu,1)[1];
    double D22u23 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[10];
    double D2Nm1u2N = D21_coeffs(RxL,0,0,0,h,D,nu,3)[3];
    
    Eigen::VectorXd d12 = D22u23*a0;
    d12(0) = D20u21;
    d12(1) = D21u22;
    d12(Ny-2) = D2Nm1u2N;
    
    D10u11 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,4)[4];
    D11u12 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,1)[1];
    D12u13 = D12_coeffs(RLy,0,0,0,h,D,nu,1)[1];
    D1Nm1u1N = D11_coeffs(RLy,RxL,0,0,h,D,nu,3)[3];
    
    
    Eigen::VectorXd d1M1 = D12u13*a1;
    d1M1(0) = D10u11;
    d1M1(1) = D11u12;
    d1M1(Ny-2) = D1Nm1u1N;
    
    D00u01 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,3)[3];
    D01u02 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,4)[4];
    D02u03 = D02_coeffs(KLy,RLy,0,0,h,D,nu,4)[4];
    D0Nm1u0N = D01_coeffs(KLy,RLy,RxL,0,h,D,nu,3)[3];
    
    
    Eigen::VectorXd d1M = D02u03*a1;
    d1M(0) = D00u01;
    d1M(1) = D01u02;
    d1M(Ny-2) = D0Nm1u0N;
    
    Eigen::VectorXd D1(Ny*Nx);
    D1 << d10,0,d11,0,d12.replicate(Nx-4,1),d1M1,0,d1M,0;
    
    //// dp2   // pad zeros at the start
    double D00u02 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,4)[4];
    double D01u03 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,5)[5];
    double D02u04 = D02_coeffs(K0y,R0y,0,0,h,D,nu,5)[5];
    
    Eigen::VectorXd d20 = D02u04*a2;
    d20(0) = D00u02;
    d20(1) = D01u03;
    
    double D10u12 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,5)[5];
    double D11u13 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,2)[2];
    double D12u14 = D12_coeffs(R0y,0,0,0,h,D,nu,2)[2];
    
    Eigen::VectorXd d21 = D12u14*a2;
    d21(0) = D10u12;
    d21(1) = D11u13;
    
    double D20u22 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,2)[2];
    double D21u23 = D21_coeffs(Rx0,0,0,0,h,D,nu,2)[2];
    double D22u24 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[12];
    
    Eigen::VectorXd d22 = D22u24*a0;
    d22(0) = D20u22;
    d22(1) = D21u23;
    
    D10u12 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,5)[5];
    D11u13 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,2)[2];
    D12u14 = D12_coeffs(RLy,0,0,0,h,D,nu,2)[2];
    
    
    Eigen::VectorXd d2M1 = D12u14*a2;
    d2M1(0) = D10u12;
    d2M1(1) = D11u13;
    
    D00u02 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,4)[4];
    D01u03 = D01_coeffs(KLy,RLy,Rx0,0,h,D,nu,5)[5];
    D02u04 = D02_coeffs(KLy,RLy,0,0,h,D,nu,5)[5];
    
    Eigen::VectorXd d2M = D02u04*a2;
    d2M(0) = D00u02;
    d2M(1) = D01u03;
    
    
    Eigen::VectorXd D2(Ny*Nx);
    D2 << d20,0,0,d21,0,0,d22.replicate(Nx-4,1),d2M1,0,0,d2M,0,0;
    //// dpNym1 // pad zeros at the start
    D01u10 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,7)[7];
    D02u11 = D02_coeffs(K0y,R0y,0,0,h,D,nu,8)[8];
    D0Nu1Nm1 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,5)[5];
    D0Nm1u1Nm2 = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,6)[6];
    
    Eigen::VectorXd dpNym10 = D02u11*a1;
    dpNym10(0) = D01u10;
    dpNym10(Ny-1-2) = D0Nm1u1Nm2;
    dpNym10(Ny-2) = D0Nu1Nm1;
    
    D11u20 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,8)[8];
    D12u21 = D12_coeffs(R0y,0,0,0,h,D,nu,9)[9];
    D1Nu2Nm1 = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,6)[6];
    D1Nm1u2Nm2 = D11_coeffs(R0y,RxL,0,0,h,D,nu,7)[7];
    
    Eigen::VectorXd dpNym11 = D12u21*a1;
    dpNym11(0) = D11u20;
    dpNym11(Ny-1-2) = D1Nm1u2Nm2;
    dpNym11(Ny-2) = D1Nu2Nm1;
    
    double D21u30 = D21_coeffs(Rx0,0,0,0,h,D,nu,9)[9];
    double D22u31 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[3];
    double D2Nu3Nm1 = D20_coeffs(KxL,RxL,0,0,h,D,nu,7)[7];
    double D2Nm1u3Nm2 = D21_coeffs(RxL,0,0,0,h,D,nu,8)[8];
    
    Eigen::VectorXd dpNym12  = D22u31*a0;
    dpNym12(1) = D21u30;
    dpNym12(Ny-2) = D2Nm1u3Nm2;
    dpNym12(Ny-1) = D2Nu3Nm1;
    
    D11u00 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,9)[9];
    D12u01 = D12_coeffs(RLy,0,0,0,h,D,nu,10)[10];
    D1Nu0Nm1 = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,7)[7];
    D1Nm1u0Nm2 = D11_coeffs(RLy,RxL,0,0,h,D,nu,10)[10];
    
    Eigen::VectorXd dpNym1M = D12u01*a1;
    dpNym1M(0) = D11u00;
    dpNym1M(Ny-1-2) = D1Nm1u0Nm2;
    dpNym1M(Ny-2) = D1Nu0Nm1;
    
    
    Eigen::VectorXd DNym1(Ny*Nx);
    DNym1 << 0,dpNym10,0,dpNym11,dpNym12.replicate(Nx-3,1),0,dpNym1M;
    //// dpNy   // pad zeros at the start
    D00u10 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,1)[1];
    D01u11 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,1)[1];
    D02u12 = D02_coeffs(K0y,R0y,0,0,h,D,nu,1)[1];
    D0Nu1N = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,1)[1];
    D0Nm1u1Nm1 = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,1)[1];
    
    Eigen::VectorXd dpNy0 = D02u12*a0;
    dpNy0(0) = D00u10;
    dpNy0(1) = D01u11;
    dpNy0(Ny-2) = D0Nm1u1Nm1;
    dpNy0(Ny-1) = D0Nu1N;
    
    D10u20 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,1)[1];
    D11u21 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,5)[5];
    D12u22 = D12_coeffs(R0y,0,0,0,h,D,nu,6)[6];
    D1Nu2N = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,1)[1];
    D1Nm1u2Nm1 = D11_coeffs(R0y,RxL,0,0,h,D,nu,5)[5];
    ////
    Eigen::VectorXd dpNy1 = D12u22*a0;
    dpNy1(0) = D10u20;
    dpNy1(1) = D11u21;
    dpNy1(Ny-2) = D1Nm1u2Nm1;
    dpNy1(Ny-1) = D1Nu2N;
    
    
    double D20u30 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,4)[4];
    double D21u31 = D21_coeffs(Rx0,0,0,0,h,D,nu,5)[5];
    double D22u32 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[7];
    double D2Nu3N = D20_coeffs(KxL,RxL,0,0,h,D,nu,4)[4];
    double D2Nm1u3Nm1 = D21_coeffs(RxL,0,0,0,h,D,nu,5)[5];
    
    Eigen::VectorXd dpNy2 = D22u32*a0;
    dpNy2(0) = D20u30;
    dpNy2(1) = D21u31;
    dpNy2(Ny-2) = D2Nm1u3Nm1;
    dpNy2(Ny-1) = D2Nu3N;
    
    D10u00 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,3)[3];
    D11u01 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,4)[4];
    D12u02 = D12_coeffs(RLy,0,0,0,h,D,nu,5)[5];
    D1Nu0N = D10_coeffs(RLy,KxL,RxL,0,h,D,nu,3)[3];
    D1Nm1u0Nm1 = D11_coeffs(RLy,RxL,0,0,h,D,nu,4)[4];
    
    Eigen::VectorXd dpNyM = D12u02*a0;
    dpNyM(0) = D10u00;
    dpNyM(1) = D11u01;
    dpNyM(Ny-2) = D1Nm1u0Nm1;
    dpNyM(Ny-1) = D1Nu0N;
    
    Eigen::VectorXd DNy(Ny*Nx);
    DNy << dpNy0,dpNy1,dpNy2.replicate((Nx-3),1),dpNyM;
    
    //// dpNyp1 // pad zeros at the start
    D00u11 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,5)[5];
    D01u12 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,6)[6];
    D02u13 = D02_coeffs(K0y,R0y,0,0,h,D,nu,7)[7];
    D0Nm1u1N = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,7)[7];
    
    Eigen::VectorXd dpNy10 = D02u13*a1 ;
    dpNy10(0) = D00u11;
    dpNy10(1) = D01u12;
    dpNy10(Ny-2) = D0Nm1u1N;
    
    D10u21 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,6)[6];
    D11u22 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,7)[7];
    D12u23 = D12_coeffs(R0y,0,0,0,h,D,nu,8)[8];
    D1Nm1u2N = D11_coeffs(R0y,RxL,0,0,h,D,nu,8)[8];
    
    
    Eigen::VectorXd dpNy11 = D12u23*a1 ;
    dpNy11(0) = D10u21;
    dpNy11(1) = D11u22;
    dpNy11(Ny-2) = D1Nm1u2N;
    
    double D20u31 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,7)[7];
    double D21u32 = D21_coeffs(Rx0,0,0,0,h,D,nu,8)[8];
    double D22u33 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[11];
    double D2Nm1u3N = D21_coeffs(RxL,0,0,0,h,D,nu,9)[9];
    
    Eigen::VectorXd dpNy12 = D22u33*a0;
    dpNy12(0) = D20u31;
    dpNy12(1) = D21u32;
    dpNy12(Ny-2) = D2Nm1u3N;
    
    D10u01 = D10_coeffs(RLy,Kx0,Rx0,0,h,D,nu,7)[7];
    D11u02 = D11_coeffs(RLy,Rx0,0,0,h,D,nu,10)[10];
    D12u03 = D12_coeffs(RLy,0,0,0,h,D,nu,11)[11];
    D1Nm1u0N = D11_coeffs(RLy,RxL,0,0,h,D,nu,9)[9];
    
    Eigen::VectorXd dpNy1M = D12u03*a1;
    dpNy1M(0) = D10u01;
    dpNy1M(1) = D11u02;
    dpNy1M(Ny-2) = D1Nm1u0N;
    
    
    Eigen::VectorXd DNy1(Ny*Nx);
    DNy1 << dpNy10,0,dpNy11,0,dpNy12.replicate(Nx-3,1),dpNy1M,0;
    //// dp2Ny  // pad zeros at the start
    D00u20 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,2)[2];
    D01u21 = D01_coeffs(K0y,R0y,Rx0,0,h,D,nu,2)[2];
    D02u22 = D02_coeffs(K0y,R0y,0,0,h,D,nu,2)[2];
    D0Nu2N = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,2)[2];
    D0Nm1u2Nm1 = D01_coeffs(K0y,R0y,RxL,0,h,D,nu,2)[2];
    
    Eigen::VectorXd dp2Ny0 = D02u22*a0;
    dp2Ny0(0) = D00u20;
    dp2Ny0(1) = D01u21;
    dp2Ny0(Ny-2) = D0Nm1u2Nm1;
    dp2Ny0(Ny-1) = D0Nu2N;
    
    D10u30 = D10_coeffs(R0y,Kx0,Rx0,0,h,D,nu,2)[2];
    D11u31 = D11_coeffs(R0y,Rx0,0,0,h,D,nu,6)[6];
    D12u32 = D12_coeffs(R0y,0,0,0,h,D,nu,7)[7];
    D1Nu3N = D10_coeffs(R0y,KxL,RxL,0,h,D,nu,2)[2];
    D1Nm1u3Nm1 = D11_coeffs(R0y,RxL,0,0,h,D,nu,6)[6];
    
    Eigen::VectorXd dp2Ny1 = D12u32*a0;
    dp2Ny1(0) = D10u30;
    dp2Ny1(1) = D11u31;
    dp2Ny1(Ny-2) = D1Nm1u3Nm1;
    dp2Ny1(Ny-1) = D1Nu3N;
    
    double D20u40 = D20_coeffs(Kx0,Rx0,0,0,h,D,nu,5)[5];
    double D21u41 = D21_coeffs(Rx0,0,0,0,h,D,nu,6)[6];
    double D22u42 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[8];
    double D2Nu4N = D20_coeffs(KxL,RxL,0,0,h,D,nu,5)[5];
    double D2Nm1u4Nm1 = D21_coeffs(RxL,0,0,0,h,D,nu,6)[6];
    
    Eigen::VectorXd dp2Ny2 = D22u42*a0;
    dp2Ny2(0) = D20u40;
    dp2Ny2(1) = D21u41;
    dp2Ny2(Ny-2) = D2Nm1u4Nm1;
    dp2Ny2(Ny-1) = D2Nu4N;
    
    Eigen::VectorXd D2Ny(Ny*Nx);
    D2Ny << dp2Ny0,dp2Ny1,dp2Ny2.replicate(Nx-3,1), Eigen::VectorXd(Ny);
    // Diag Biharmonic
    
    std::vector<Eigen::VectorXd> BHdiags = {
        Dm2Ny,
        DmNym1, DmNy, DmNy1,
        Dm2, Dm1, D0, D1, D2,
        DNym1, DNy, DNy1,
        D2Ny
    };
    
    std::vector<int>dn = {-(2*(Ny)),-(Ny+1),-(Ny),-(Ny-1), -2,1,0,1,2, (Ny-1),(Ny),(Ny+1), 2*(Ny)};
    
    spdiags(BHdiags, dn, biharm);
    double oh = 1/h;
    biharm *= oh*oh*oh*oh;
}

#endif /* bhmat_h */
