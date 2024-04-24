//
//  spdiags.h
//  using_eigen
//
//  Created by Matthew Hamilton on 19/04/2024.
//

#ifndef spdiags_h
#define spdiags_h

#include "range.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include "Eigen/Eigen"

using util::lang::indices;

void spdiags(Eigen::SparseMatrix<double> A)
{
    
//    //Extraction of the diagnols before the main diagonal
//    std::vector<double> vec1; 
//    int flag=0;
//    int l=0;
//    int j=0;
//    std::vector<std::vector<double> > diagD;
//    std::vector<std::vector<double> > diagG;  
//    int z=0;
//    
//    
//    for(int i=0;i<A.rows();i++)
//    {
//        int l=i;
//        do
//        {
//            if(A.coeff(l,j)!=0)
//                flag=1;
//            
//            vec1.push_back(A.coeff(l,j));
//            l++;j++;
//            
//        }while(l<A.rows() && j<A.cols());
//        
//        if(flag==1) {diagG.resize(diagG.size()+1);diagG[z]=vec1; z++; }
//        vec1.clear(); l=0;j=0; flag=0; std::cout<<std::endl;
//    }
//    
//    
//    flag=0;z=0; vec1.clear();
//    
//    // Extraction of the diagonals after the main diagonal
//    for(int i=1;i<A.cols();i++)
//    {
//        int l=i;
//        do
//        {
//            if(A.coeff(j,l)!=0)
//                flag=1;
//            
//            vec1.push_back(A.coeff(j,l));
//            l++;j++;
//            
//        }while(l<A.cols() && j<A.rows());
//        
//        if(flag==1) {diagD.resize(diagD.size()+1);diagD[z]=vec1; z++;  }
//        vec1.clear(); l=0;j=0; flag=0; std::cout<<std::endl;
//    }
//    // End extraction of the diagonals
//    
//    Eigen::VectorXi d = Eigen::VectorXi::Zero(A.rows() + A.cols() - 1);
//    
//    for (int k=0; k < A.outerSize(); ++k)
//    {
//        for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
//        {
//            d(it.col() - it.row() + A.rows() - 1) = 1;
//            
//        }
//    }
//    
//    
//    int num_diags = d.sum();
//    Eigen::MatrixXd B(std::min(A.cols(), A.rows()), num_diags);
//    
//    // fill B with diagonals
//    Eigen::ArrayXd v;
//    int B_col_idx = 0;
//    int B_row_sign = A.rows() >= A.cols() ? 1 : -1;
//    int indG=diagG.size()-1; int indD=0;
//    
//    for (int i = 1 - A.rows(); i <=A.cols() - 1; i++)
//    {
//        if (d(i + A.rows() - 1))
//        {
//            if(i<1)
//            {   v.resize(diagG[indG].size());
//                for(int i=0;i<diagG[indG].size();i++)
//                {
//                    v(i)=diagG[indG][i];
//                }
//                
//                int B_row_start = std::max(0, B_row_sign * i);
//                B.block(B_row_start, B_col_idx, diagG[indG].size(), 1) = v;
//                
//                B_col_idx++;
//                indG--;
//            }
//            else
//            {
//                v.resize(diagD[indD].size());
//                
//                for(int i=0;i<diagD[indD].size();i++)
//                {
//                    v(i)=diagD[indD][i] ;
//                }
//                
//                int B_row_start = std::max(0, B_row_sign * i);
//                B.block(B_row_start, B_col_idx, diagD[indD].size(), 1) = v;
//                
//                B_col_idx++;
//                indD++;
//            }
//        }
//        
//    }
//    std::cout<<B<<std::endl; //the result of the function
}//end of the function



extern void vector_list(std::vector<int>list)
{
    for (auto a: list) {
        std::cout << a << std::endl;
    }
}

void spdiags(std::vector<Eigen::VectorXd> diagVectors,
             std::vector<int> diagNums,
             Eigen::SparseMatrix<double>& A)
{
    std::vector<Eigen::Triplet<double>> coefficients;
                    
    for(auto i : indices(diagVectors))
    {
        auto& diag = diagVectors.at(i);
        auto& diagNum = diagNums.at(i);
        
        if (diagNum == 0) {
            for (int i = 0; i < diag.size(); i ++)
                coefficients.push_back(Eigen::Triplet<double>(i,i,diag[i]));
        }
        else if (diagNum > 0) {
            for (int i = 0; i < (diag.size()-diagNum); i ++)
                coefficients.push_back(Eigen::Triplet<double>(i,i+diagNum,diag[i]));
        }
        else
        {            
            for (int i = 0; i < (diag.size()-std::abs(diagNum)); i ++)
                coefficients.push_back(Eigen::Triplet<double>(i+std::abs(diagNum),i,diag[i]));
        }
    }

    A.setFromTriplets(coefficients.begin(), coefficients.end());
}

#endif /* spdiags_h */


