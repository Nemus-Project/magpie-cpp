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
    int Nx = 5;
    int Ny = 5;

    Eigen::SparseMatrix<double> d(Ny*Nx,Ny*Nx);
    bhmat(d,{0,0,0,0,0,0,0,0}, {Nx,Ny}, 0.1, 0.1, 12, 0.3);
    
    return 0;
}
