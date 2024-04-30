//
//  main.cpp
//  using_eigen
//
//  Created by Matthew Hamilton on 18/04/2024.
//

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

#include "range.h"
#include "spdiags.h"
#include "D_Coeffs.h"
#include "bhmat.h"

int main(int argc, const char * argv[]) 
{
    int Nx = 20;
    int Ny = 20;

    Eigen::SparseMatrix<double> biharm(Ny*Nx,Ny*Nx);
    bhmat(biharm,{0,0,0,0,0,0,0,0}, {Nx,Ny}, 0.01, 0.01, 200e9, 0.3);
    
    Eigen::EigenSolver<Eigen::SparseMatrix<double>> eigensolver;
//    eigensolver.compute(biharm);
    
//    eigensolver.info();
    
//    switch (info) {
//        case Eigen::Success:
//            std::cout << "Success" << std::endl;
//            break;
//        case Eigen::NumericalIssue:
//            std::cout << "NumericalIssue" << std::endl;
//            break;
//        case Eigen::NoConvergence:
//            std::cout << "NoConvergence" << std::endl;
//            break;
//        case Eigen::InvalidInput:
//            std::cout << "InvalidInput" << std::endl;
//            break;
//    }
////    solver.factorize(biharm);
////    solver.compute(biharm);
////    Eigen::VectorXd x = solver.solve(b);
//    
//    
//    std::cout << x << std::endl;
    
    
    return 0;
}
