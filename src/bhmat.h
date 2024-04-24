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

using util::lang::indices;

void bhmat(std::array<double,8>BCs, std::array<int,2> Nxy, double h, double Lz, double E, double nu)
{
    const auto&  [K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = BCs;
    const int Nx = Nxy[0];
    const int Ny = Nxy[1];
    const double D = E * (Lz*Lz*Lz) / 12 / (1-(nu*nu));
    
    Eigen::SparseMatrix<double> biharm(Nx,Ny);
    Eigen::VectorXd a0(Ny+1);
    Eigen::VectorXd a1(Ny);
    Eigen::VectorXd a2(Ny-1);
    
    //// dm2Ny  // pad zeros at the end
    double D20u00 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,6)[6];
    double D21u01 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,7)[7];
    double D22u02 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,4)[4];
    double D2Nu0N = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,6)[6];
    double D2Nm1u0Nm1 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,7)[7];
    
    //(d\w+) = (D\w+)\*a(\d);
    //(d\w+) = (D\w+)*a(\d);
    /*
     for(auto i: indices(std::vector<int>{$2}))
     {
     $1(std::vector<int>{$2}[i]) = std::vector<int>{$3}[i];
     }
     */
    Eigen::VectorXd  dm2Ny0 = D22u02*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dm2Ny0(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D10u30 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,2)[2];
    double D11u31 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,6)[6];
    double D12u32 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,7)[7];
    double D1Nu3N = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,2)[2];
    double D1Nm1u3Nm1 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,6)[6];
    
    Eigen::VectorXd dm2Ny1 = D12u32*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dm2Ny1(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D00u20 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,2)[2];
    double D01u21 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,2)[2];
    double D02u22 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,2)[2];
    double D0Nu2N = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,2)[2];
    double D0Nm1u2Nm1 = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,2)[2];
    
    Eigen::VectorXd dm2Ny2 = D02u22*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dm2Ny2(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd Dm2Ny(Ny*Nx);
    
    Dm2Ny << dm2Ny0.replicate(Nx-3,1), dm2Ny1, dm2Ny2;
    
    //// dmNym1 // pad zeros at the end
    double D11u00 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,9)[9];
    double D12u01 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,10)[10];
    double D1Nu0Nm1 = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,7)[7];
    double D1Nm1u0Nm2 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,10)[10];
    
    
    Eigen::VectorXd dmNym10 = D12u01*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dmNym10(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D21u10 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,10)[10];
    double D22u11 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[1];
    double D2Nu1Nm1 = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,8)[8];
    double D2Nm1u1Nm2 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,11)[11];
    
    Eigen::VectorXd dmNym11 = D22u11*a0;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dmNym11(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D11u20 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,8)[8];
    double D12u21 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,9)[9];
    double D1Nm1u2Nm2 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,7)[7];
    double D1Nu2Nm1 = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,6)[6];
    
    Eigen::VectorXd dmNym1M1 = D12u21*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dmNym1M1(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D01u10 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,7)[7];
    double D02u11 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,8)[8];
    double D0Nu1Nm1 = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,5)[5];
    double D0Nm1u1Nm2 = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,6)[6];
    
    Eigen::VectorXd dmNym1M = D02u11*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dmNym1M(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    //    VectorXd vec_joined(vec1.size() + vec2.size());
    //    vec_joined << vec1, vec2;
    Eigen::VectorXd DmNym1(Ny*Nx);
    DmNym1 << dmNym10,0,dmNym11.replicate(Nx-3,1),dmNym1M1,0,dmNym1M;
    //    DmNym1 = [DmNym1; zeros((Ny + 2),1)];
    //// dmNy   // pad zeros at the end
    
    double D10u00 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,3)[3];
    double D11u01 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,4)[4];
    double D12u02 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,5)[5];
    double D1Nu0N = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,3)[3];
    double D1Nm1u0Nm1 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,4)[4];
    
    Eigen::VectorXd dmNy0 = D12u02*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dmNy0(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D20u10 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,3)[3];
    double D21u11 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,4)[4];
    double D22u12 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,5)[5];
    double D2Nu1N = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,3)[3];
    double D2Nm1u1Nm1 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,4)[4];
    
    Eigen::VectorXd dmNy1 = D22u12*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dmNy1(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D10u20 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,1)[1];
    double D11u21 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,5)[5];
    double D12u22 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,6)[6];
    double D1Nu2N = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,1)[1];
    double D1Nm1u2Nm1 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,5)[5];
    
    Eigen::VectorXd dmNyM1 = D12u22*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dmNyM1(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D00u10 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,1)[1];
    double D01u11 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,1)[1];
    double D02u12 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,1)[1];
    double D0Nu1N = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,1)[1];
    double D0Nm1u1Nm1 = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,1)[1];
    
    Eigen::VectorXd dmNyM = D02u12*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dmNyM(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd DmNy(Ny*Nx);
    DmNy << dmNy0,dmNy1.replicate((Nx-3),1),dmNyM1,dmNyM;
    //    DmNy = [DmNy; zeros((Ny + 1),1)];
    ////assert(all(abs(diag(biHarm,Ny+1) - DmNy) <= eps), "D0 incorrect");
    
    
    //// dmNyp1 // pad zeros at the
    double D10u01 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,7)[7];
    double D11u02 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,10)[10];
    double D12u03 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,11)[11];
    double D1Nm1u0N = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,9)[9];
    
    Eigen::VectorXd dmNy10 = D12u03*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dmNy10(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D20u11 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,8)[8];
    double D21u12 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,11)[11];
    double D22u13 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,9)[9];
    double D2Nm1u1N = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,10)[10];
    
    
    Eigen::VectorXd dmNy11 = D22u13*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dmNy11(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D10u21 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,6)[6];
    double D11u22 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,7)[7];
    double D12u23 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,8)[8];
    double D1Nm1u2N = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,8)[8];
    
    
    Eigen::VectorXd dmNy1M1 = D12u23*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dmNy1M1(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D00u11 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,5)[5];
    double D01u12 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,6)[6];
    double D02u13 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,7)[7];
    double D0Nm1u1N = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,7)[7];
    
    Eigen::VectorXd dmNy1M = D02u13*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dmNy1M(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    
    Eigen::VectorXd DmNy1(Ny*Nx);
    
    DmNy1 << 0,dmNy10,dmNy11.replicate(Nx-3,1),0,dmNy1M1,0,dmNy1M,0;
    //    DmNy1 = [DmNy1; zeros((Ny),1) ];
    //// dm2   // pad zeros at the end
    double D02u00 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,6)[6];
    double D0Nu0Nm2 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,4)[4];
    double D0Nm1u0Nm3 = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,5)[5];
    
    
    Eigen::VectorXd dm20 = D02u00*a2;
    for(auto i: indices(std::vector<int>{Ny-2,Ny-1}))
    {
        dm20(std::vector<int>{Ny-2,Ny-1}[i]) = std::vector<double>{}[i];
    }
    
    double D12u10 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,4)[4];
    double D1Nu1Nm2 = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,5)[5];
    double D1Nm1u1Nm3 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,2)[2];
    
    Eigen::VectorXd dm21 = D12u10*a2;
    for(auto i: indices(std::vector<int>{Ny-2,Ny-1}))
    {
        dm21(std::vector<int>{Ny-2,Ny-1}[i]) = std::vector<double>{}[i];
    }
    
    double D22u20 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,0)[0];
    double D2Nu2Nm2 = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,2)[2];
    double D2Nm1u2Nm3 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,2)[2];
    
    
    Eigen::VectorXd dm22 = D22u20*a0;
    for(auto i: indices(std::vector<int>{Ny-2,Ny-1}))
    {
        dm22(std::vector<int>{Ny-2,Ny-1}[i]) = std::vector<double>{}[i];
    }
    
    D12u10 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,4)[4];
    D1Nu1Nm2 = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,5)[5];
    D1Nm1u1Nm3 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,2)[2];
    
    
    Eigen::VectorXd dm2M1 = D12u10*a2 ;
    for(auto i: indices(std::vector<int>{Ny-2,Ny-1}))
    {
        dm2M1(std::vector<int>{Ny-2,Ny-1}[i]) = std::vector<double>{}[i];
    }
    
    D02u00 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,6)[6];
    D0Nu0Nm2 = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,4)[4];
    D0Nm1u0Nm3 = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,5)[5];
    
    Eigen::VectorXd dm2M = D02u00*a2;
    for(auto i: indices(std::vector<int>{Ny-2,Ny-1}))
    {
        dm2M(std::vector<int>{Ny-2,Ny-1}[i]) = std::vector<double>{}[i];
    }
    
    
    Eigen::VectorXd Dm2(Ny*Nx);
    Dm2 << dm20,0,0,dm21,0,0,dm22.replicate(Nx+1-4,1),dm2M1,0,0,dm2M;
    //    Dm2 = [Dm2,0,0];
    //// dm1   // pad zeros at the end
    
    
    double D01u00 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,3)[3];
    double D02u01 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,3)[3];
    double D0Nu0Nm1 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,3)[3];
    double D0Nm1u0Nm2 = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,4)[4];
    
    Eigen::VectorXd dm10 = D02u01*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dm10(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D11u10 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,3)[3];
    double D12u11 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,3)[3];
    double D1Nu1Nm1 = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,4)[4];
    double D1Nm1u1Nm2 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,1)[1];
    
    
    Eigen::VectorXd dm11 = D12u11*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dm11(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D21u20 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,3)[3];
    double D22u21 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,2)[2];
    double D2Nu2Nm1 = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,1)[1];
    double D2Nm1u2Nm2 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,1)[1];
    
    
    Eigen::VectorXd dm12 = D22u21*a0;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dm12(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D11u10 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,3)[3];
    D12u11 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,3)[3];
    D1Nm1u1Nm2 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,1)[1];
    D1Nu1Nm1 = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,4)[4];
    
    
    Eigen::VectorXd dm1M1 = D12u11*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dm1M1(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D01u00 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,3)[3];
    D02u01 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,3)[3];
    D0Nm1u0Nm2 = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,4)[4];
    D0Nu0Nm1 = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,3)[3];
    
    
    Eigen::VectorXd dm1M = D02u01*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dm1M(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd Dm1(Ny*Nx);
    Dm1 << dm10,0,dm11,0,dm12.replicate(Nx+1-4,1),dm1M1,0,dm1M;
    //    Dm1 = [Dm1,0];
    //// d00
    
    double D00u00 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,0)[0];
    double D01u01 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,0)[0];
    double D02u02 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,0)[0];
    double D0Nu0N = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,0)[0];
    double D0Nm1u0Nm1 = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,0)[0];
    
    Eigen::VectorXd d00 = D02u02*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        d00(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D10u10 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,0)[0];
    double D11u11 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,0)[0];
    double D12u12 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,0)[0];
    double D1Nu1N = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,0)[0];
    double D1Nm1u1Nm1 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d01 = D12u12*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        d01(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D20u20 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,0)[0];
    double D21u21 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,0)[0];
    double D22u22 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,6)[6];
    double D2Nu2N = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,0)[0];
    double D2Nm1u2Nm1 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d02 = D22u22*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        d02(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    
    D10u10 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,0)[0];
    D11u11 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,0)[0];
    D12u12 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,0)[0];
    D1Nm1u1Nm1 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,0)[0];
    D1Nu1N = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d0Mm = D12u12*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        d0Mm(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    D00u00 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,0)[0];
    D01u01 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,0)[0];
    D02u02 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,0)[0];
    D0Nm1u0Nm1 = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,0)[0];
    D0Nu0N = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,0)[0];
    
    
    Eigen::VectorXd d0M = D02u02*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        d0M(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd D0(Ny*Nx);
    D0 << d00,d01,d02.replicate((Nx+1-4),1),d0Mm,d0M;
    
    // //assert(all(abs(diag(biHarm,0) - D0) <= eps), "D0 incorrect");
    
    //// dp1   // pad zeros at the start
    double D00u01 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,3)[3];
    double D01u02 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,4)[4];
    double D02u03 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,4)[4];
    double D0Nm1u0N = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,3)[3];
    
    Eigen::VectorXd d10 = D02u03*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        d10(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D10u11 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,4)[4];
    double D11u12 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,1)[1];
    double D12u13 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,1)[1];
    double D1Nm1u1N = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,3)[3];
    
    Eigen::VectorXd d11 = D12u13*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        d11(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D20u21 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,1)[1];
    double D21u22 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,1)[1];
    double D22u23 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,10)[10];
    double D2Nm1u2N = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,3)[3];
    
    Eigen::VectorXd d12 = D22u23*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        d12(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D10u11 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,4)[4];
    D11u12 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,1)[1];
    D12u13 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,1)[1];
    D1Nm1u1N = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,3)[3];
    
    
    Eigen::VectorXd d1M1 = D12u13*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        d1M1(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D00u01 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,3)[3];
    D01u02 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,4)[4];
    D02u03 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,4)[4];
    D0Nm1u0N = D01_coeffs(KLy,RLy,RxL,Rx0,h,D,nu,3)[3];
    
    
    Eigen::VectorXd d1M = D02u03*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        d1M(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd D1(Ny*Nx);
    D1 << d10,0,d11,0,d12.replicate(Nx+1-4,1),d1M1,0,d1M;
    
    //    D1 = [0;D1];
    //// dp2   // pad zeros at the start
    double D00u02 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,4)[4];
    double D01u03 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,5)[5];
    double D02u04 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,5)[5];
    
    Eigen::VectorXd d20 = D02u04*a2;
    for(auto i: indices(std::vector<int>{1,2}))
    {
        d20(std::vector<int>{1,2}[i]) = std::vector<double>{}[i];
    }
    
    double D10u12 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,5)[5];
    double D11u13 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,2)[2];
    double D12u14 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,2)[2];
    
    Eigen::VectorXd d21 = D12u14*a2;
    for(auto i: indices(std::vector<int>{1,2}))
    {
        d21(std::vector<int>{1,2}[i]) = std::vector<double>{}[i];
    }
    
    double D20u22 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,2)[2];
    double D21u23 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,2)[2];
    double D22u24 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,12)[12];
    
    Eigen::VectorXd d22 = D22u24*a0;
    for(auto i: indices(std::vector<int>{1,2}))
    {
        d22(std::vector<int>{1,2}[i]) = std::vector<double>{}[i];
    }
    
    D10u12 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,5)[5];
    D11u13 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,2)[2];
    D12u14 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,2)[2];
    
    Eigen::VectorXd d2M1 = D12u14*a2;
    for(auto i: indices(std::vector<int>{1,2}))
    {
        d2M1(std::vector<int>{1,2}[i]) = std::vector<double>{}[i];
    }
    
    D00u02 = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,4)[4];
    D01u03 = D01_coeffs(KLy,RLy,Rx0,Rx0,h,D,nu,5)[5];
    D02u04 = D02_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,5)[5];
    
    Eigen::VectorXd d2M = D02u04*a2;
    for(auto i: indices(std::vector<int>{1,2}))
    {
        d2M(std::vector<int>{1,2}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd D2(Ny*Nx);
    D2 << d20,0,0,d21,0,0,d22.replicate(Nx+1-4,1),d2M1,0,0,d2M;
    
    //    D2 = [0,0,D2];
    //// dpNym1 // pad zeros at the start
    D01u10 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,7)[7];
    D02u11 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,8)[8];
    D0Nu1Nm1 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,5)[5];
    D0Nm1u1Nm2 = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,6)[6];
    
    Eigen::VectorXd dpNym10 = D02u11*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dpNym10(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D11u20 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,8)[8];
    D12u21 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,9)[9];
    D1Nu2Nm1 = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,6)[6];
    D1Nm1u2Nm2 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,7)[7];
    
    Eigen::VectorXd dpNym11 = D12u21*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dpNym11(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D21u30 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,9)[9];
    double D22u31 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,3)[3];
    double D2Nu3Nm1 = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,7)[7];
    double D2Nm1u3Nm2 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,8)[8];
    
    Eigen::VectorXd dpNym12  = D22u31*a0;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dpNym12(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D11u00 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,9)[9];
    D12u01 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,10)[10];
    D1Nm1u0Nm2 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,10)[10];
    D1Nu0Nm1 = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,7)[7];
    
    Eigen::VectorXd dpNym1M = D12u01*a1;
    for(auto i: indices(std::vector<int>{1,Ny-1,Ny}))
    {
        dpNym1M(std::vector<int>{1,Ny-1,Ny}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd DNym1(Ny*Nx);
    DNym1 << 0,dpNym10,0,dpNym11,dpNym12.replicate(Nx-3,1),0,dpNym1M,0;
    //    DNym1 = [zeros((Ny),1);DNym1];
    //// dpNy   // pad zeros at the start
    D00u10 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,1)[1];
    D01u11 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,1)[1];
    D02u12 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,1)[1];
    D0Nu1N = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,1)[1];
    D0Nm1u1Nm1 = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,1)[1];
    
    Eigen::VectorXd dpNy0 = D02u12*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dpNy0(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    D10u20 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,1)[1];
    D11u21 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,5)[5];
    D12u22 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,6)[6];
    D1Nu2N = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,1)[1];
    D1Nm1u2Nm1 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,5)[5];
    ///
    Eigen::VectorXd dpNy1 = D12u22*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dpNy1(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    
    double D20u30 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,4)[4];
    double D21u31 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,5)[5];
    double D22u32 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,7)[7];
    double D2Nu3N = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,4)[4];
    double D2Nm1u3Nm1 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,5)[5];
    
    Eigen::VectorXd dpNy2 = D22u32*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dpNy2(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    D10u00 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,3)[3];
    D11u01 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,4)[4];
    D12u02 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,5)[5];
    D1Nu0N = D10_coeffs(RLy,KxL,RxL,Rx0,h,D,nu,3)[3];
    D1Nm1u0Nm1 = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,4)[4];
    
    Eigen::VectorXd dpNyM = D12u02*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dpNyM(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    Eigen::VectorXd DNy(Ny*Nx);
    DNy<<dpNy0,dpNy1,dpNy2.replicate((Nx-3),1),dpNyM;
//    DNy = [zeros((Ny+1),1);DNy];
    
    //// dpNyp1 // pad zeros at the start
    D00u11 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,5)[5];
    D01u12 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,6)[6];
    D02u13 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,7)[7];
    D0Nm1u1N = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,7)[7];
    
    Eigen::VectorXd dpNy10 = D02u13*a1 ;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dpNy10(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D10u21 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,6)[6];
    D11u22 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,7)[7];
    D12u23 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,8)[8];
    D1Nm1u2N = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,8)[8];
    
    
    Eigen::VectorXd dpNy11 = D12u23*a1 ;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dpNy11(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    double D20u31 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,7)[7];
    double D21u32 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,8)[8];
    double D22u33 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,11)[11];
    double D2Nm1u3N = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,9)[9];
    
    Eigen::VectorXd dpNy12 = D22u33*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dpNy12(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    D10u01 = D10_coeffs(RLy,Kx0,Rx0,Rx0,h,D,nu,7)[7];
    D11u02 = D11_coeffs(RLy,Rx0,Kx0,Rx0,h,D,nu,10)[10];
    D12u03 = D12_coeffs(RLy,R0y,Kx0,Rx0,h,D,nu,11)[11];
    D1Nm1u0N = D11_coeffs(RLy,RxL,Kx0,Rx0,h,D,nu,9)[9];
    
    Eigen::VectorXd dpNy1M = D12u03*a1;
    for(auto i: indices(std::vector<int>{1,2,Ny}))
    {
        dpNy1M(std::vector<int>{1,2,Ny}[i]) = std::vector<double>{}[i];
    }
    
    
    Eigen::VectorXd DNy1(Ny*Nx);
    DNy1<< dpNy10,0,dpNy11,0,dpNy12.replicate(Nx-3,1),dpNy1M;
//    DNy1 = [zeros((Ny+2),1);DNy1];
    //// dp2Ny  // pad zeros at the start
    D00u20 = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,2)[2];
    D01u21 = D01_coeffs(K0y,R0y,Rx0,Rx0,h,D,nu,2)[2];
    D02u22 = D02_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,2)[2];
    D0Nu2N = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,2)[2];
    D0Nm1u2Nm1 = D01_coeffs(K0y,R0y,RxL,Rx0,h,D,nu,2)[2];
    
    Eigen::VectorXd dp2Ny0 = D02u22*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dp2Ny0(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    D10u30 = D10_coeffs(R0y,Kx0,Rx0,Rx0,h,D,nu,2)[2];
    D11u31 = D11_coeffs(R0y,Rx0,Kx0,Rx0,h,D,nu,6)[6];
    D12u32 = D12_coeffs(R0y,R0y,Kx0,Rx0,h,D,nu,7)[7];
    D1Nu3N = D10_coeffs(R0y,KxL,RxL,Rx0,h,D,nu,2)[2];
    D1Nm1u3Nm1 = D11_coeffs(R0y,RxL,Kx0,Rx0,h,D,nu,6)[6];
    
    Eigen::VectorXd dp2Ny1 = D12u32*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dp2Ny1(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    
    double D20u40 = D20_coeffs(Kx0,Rx0,Kx0,Rx0,h,D,nu,5)[5];
    double D21u41 = D21_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,6)[6];
    double D22u42 = D22_coeffs(Rx0,R0y,Kx0,Rx0,h,D,nu,8)[8];
    double D2Nu4N = D20_coeffs(KxL,RxL,Kx0,Rx0,h,D,nu,5)[5];
    double D2Nm1u4Nm1 = D21_coeffs(RxL,R0y,Kx0,Rx0,h,D,nu,6)[6];
    
    Eigen::VectorXd dp2Ny2 = D22u42*a0;
    for(auto i: indices(std::vector<int>{1,2,Ny,Ny+1}))
    {
        dp2Ny2(std::vector<int>{1,2,Ny,Ny+1}[i]) = std::vector<double>{}[i];
    }
    Eigen::VectorXd D2Ny(Ny*Nx);
    D2Ny << dp2Ny0,dp2Ny1,dp2Ny2.replicate(Nx-3,1);
//    D2Ny = [zeros((Ny+1)*2,1);D2Ny];
    
    std::vector<Eigen::VectorXd> BHdiags = {
        Dm2Ny,
        DmNym1, DmNy, DmNy1,
        Dm2, Dm1, D0, D1, D2,
        DNym1, DNy, DNy1,
        D2Ny
    };
    
    std::vector<int>dn = {-(2*(Ny+1)),-(Ny+2),-(Ny+1),-(Ny), -2,1,0,1,2, (Ny),(Ny+1),(Ny+2), 2*(Ny+1)};
    
    spdiags(BHdiags, dn, biharm);
    double oh = 1/h;
    biharm *= oh*oh*oh*oh;
    
}

#endif /* bhmat_h */
