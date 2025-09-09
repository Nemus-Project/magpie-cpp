//
//  magpie.h
//  magpie
//
//  Created by admin on 09/09/2025.
//

#ifndef magpie_h
#define magpie_h

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
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

#include "bhmat.h"



/// <#Description#>
/// - Parameters:
///   - rho: <#rho description#>
///   - E: <#E description#>
///   - nu: <#nu description#>
///   - ldim: <#ldim description#>
///   - h: <#h description#>
///   - BCs: <#BCs description#>
void magpie(double rho, 
            double E,
            double nu, 
            std::array<int,3> ldim, 
            double h,
            std::array<double,8> BCs,
            int Nm = 0)
{
    const auto& [Lx, Ly, Lz] = ldim;
    double D = E * (Lz * Lz * Lz) / 12.0 / (1.0 - (nu * nu));

    int Nx = Lx / h;
    int Ny = Ly / h;
    if(Nm == 0)
        Nm = 16;
    
    
    Eigen::SparseMatrix<double> biharm(Ny*Nx,Ny*Nx);
    bhmat(biharm, BCs, {Nx, Ny}, h, Lz, E, nu);
    
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
        auto Dm = eigs.eigenvalues();
        auto Q = eigs.eigenvectors();
        std::vector<double> Om;
        Om.reserve(Nm);
        
        for (int i = 0; i < nconv; i++) {
            Om[i] = (std::sqrt(std::abs(Dm[i].real()))*std::sqrt(D/rho/Lz)) / M_2_PI;
        }
    }
}
#endif /* magpie_h */
