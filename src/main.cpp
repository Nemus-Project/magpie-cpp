//
//  main.cpp
//  using_eigen
//
//  Created by Matthew Hamilton on 18/04/2024.
//

#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <array>
#include <Eigen/Eigen>

#include "range.h"
#include "spdiags.h"
#include "D_Coeffs.h"

int main(int argc, const char * argv[]) {
    //    Eigen::MatrixXd m(2,2);
    //    m(0,0) = 3;
    //    m(1,0) = 2.5;
    //    m(0,1) = -1;
    //    m(1,1) = m(1,0) + m(0,1);
    //    std::cout << m << std::endl;
    //
    //    Eigen::MatrixXd a(2,2);
    //        a << 1,2,3,4;
    
    Eigen::Matrix< double, 3, 1> v ;
    v << 1, 2, 3;
    
    Eigen::SparseMatrix<double> b(3,3);
    
    std::cout << Eigen::MatrixXd(b) << "\n";
           
    Eigen::SparseMatrix<double> S(3, 3);
    
    std::cout << "S =\n" << S << "\n";
    
    int vecSize = 4;
    Eigen::VectorXd s(vecSize);

    std::cout <<  s << "\n";

    std::vector<Eigen::Triplet<double>> coefficients;
    
    coefficients.push_back(Eigen::Triplet<double>(1,2,4.3));
    
    std::cout << coefficients[0].value() << ' '
    << coefficients[0].col() << ' ' << coefficients[0].row()   << std::endl;
    
    v[0];
    
    Eigen::SparseMatrix<double> A(6,6);

    Eigen::VectorXd a(6);
    a.setOnes();
    
    spdiags({a}, {0}, A);

    std::cout <<  A.toDense() << "\n";
    
    return 0;
}
