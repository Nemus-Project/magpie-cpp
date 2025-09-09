//
//  main.cpp
//  using_eigen
//
//  Created by Matthew Hamilton on 18/04/2024.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/GenEigsSolver.h>
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>

#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseGenComplexShiftSolve.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseHermMatProd.h>
#include <Spectra/MatOp/SparseRegularInverse.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/SymShiftInvert.h>

#include <complex>

#include "range.h"
#include "spdiags.h"
#include "bhmat.h"
#include <iomanip>






using Complex = std::complex<double>;

int main(int argc, const char * argv[])
{
    double Lx = 1.10, Ly = 0.8, Lz = 5e-3;
    double ldim[3] = {Lx, Ly, Lz};       // plate dimensions [x, y, z] in metres
    double E    = 9.0e+9 ;        // Young's mod [Pa]
    double rho  = 8765 ;            // density [kg/m^3]
    double nu   = 0.3 ;             // poisson's ratio
    int Nm   = 20;               // number of modes to compute
    double h    = std::sqrt(Lx*Ly)*0.01; // Grid Spacing
    double BC   = 1e15;   // elastic constants around the edges
    
    double D    = E * Lz*Lz*Lz / 12.0 / (1.0-(nu*nu));
        
    int Nx = 110;
    int Ny = 80;

    Eigen::SparseMatrix<double> biharm(Ny*Nx,Ny*Nx);
    
    bhmat(biharm,{BC,BC,BC,BC,BC,BC,BC,BC}, {Nx,Ny}, h, Lz, E, nu);
    
    std::cout << "Nx: " << Nx << '\n';
    std::cout << "Ny: " << Ny << '\n';
    

    Eigen::SparseMatrix<double> I(Nx*Ny, Nx*Ny);
    I.setIdentity();
    
    Spectra::SparseCholesky<double> opB(I);

    Spectra::SparseGenRealShiftSolve<double> op(biharm);
    Spectra::GenEigsRealShiftSolver<Spectra::SparseGenRealShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0000);

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
                
        for (auto val: evalues) {
            std::cout << std::sqrt(std::abs(val.real()))*std::sqrt(D/rho/Lz) << std::endl;
        }
        
    }
    
    return 0;
}
