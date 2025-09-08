//
//  main.cpp
//  using_eigen
//
//  Created by Matthew Hamilton on 18/04/2024.
//

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include "range.h"
#include "spdiags.h"
#include "D_Coeffs.h"
#include "bhmat.h"
#include "Dxx_coeffs.h"

int main(int argc, const char * argv[]) 
{
    int Nx = 20;
    int Ny = 20;

    Eigen::SparseMatrix<double> biharm(Ny*Nx,Ny*Nx);
    
    
    double Lx = 0.10, Ly = 0.08, Lz = 0.81e-3;
    double ldim[3] = {Lx, Ly, Lz};       // plate dimensions [x, y, z] in metres
    double E    = 1.01e+11 ;        // Young's mod [Pa]
    double rho  = 8765 ;            // density [kg/m^3]
    double nu   = 0.3 ;             // poisson's ratio
    double Nm   = 16;               // number of modes to compute
    double h    = std::sqrt(Lx*Ly)*0.01; // Grid Spacing
    double BC   = 1e15;   // elastic constants around the edges
    double D    = E * Lz*Lz*Lz / 12.0 / (1.0-(nu*nu));
    
    auto K0y = BC;
    auto Kx0 = BC;
    auto KLy = BC;
    auto KxL = BC;
    auto R0y = BC;
    auto Rx0 = BC;
    auto RLy = BC;
    auto RxL = BC;
    
    auto ans1  = D01_coeffs(KLy, RLy, Rx0, h, D, nu);
//    auto ans2  = D01_coeffs_x(KLy,RLy,Kx0,Rx0,h,D,nu);
    //(BCs,Nxy,h,Lz,E,nu)
    bhmat(biharm,{BC,BC,BC,BC,BC,BC,BC,BC}, {Nx,Ny}, h, Lz, E, nu);
    
    std::cout << biharm.coeff(0,0) << std::endl;
    
    
    // Construct matrix operation object using the wrapper class SparseGenMatProd
    Spectra::SparseGenMatProd<double> op(biharm);
 
    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double>> eigs(op, 3, 6);
 
    // Initialize and compute
    eigs.init();
    auto nconv = eigs.compute(Spectra::SortRule::LargestMagn);
 

    switch (eigs.info()) {
        case Spectra::CompInfo::Successful:
            std::cout << "Success" << std::endl;
            break;
        case Spectra::CompInfo::NumericalIssue:
            std::cout << "NumericalIssue" << std::endl;
            break;
        case Spectra::CompInfo::NotConverging:
            std::cout << "NotConverging" << std::endl;
            break;
        case Spectra::CompInfo::NotComputed:
            std::cout << "NotComputed" << std::endl;
            break;
    }

    if(eigs.info() == Spectra::CompInfo::Successful) // Retrieve results
    {
        auto evalues = eigs.eigenvalues();
        auto evectors = eigs.eigenvectors();

        std::cout << "Number Converged eigenvalues:" << nconv << std::endl;
        std::cout << "Converged eigenvalues:\n" << evalues << std::endl;
        std::cout << "First eigenvector:\n" << evectors.col(0).transpose() << std::endl;
    }
    
    return 0;
}
