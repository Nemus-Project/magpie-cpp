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

#include "magpie.h"

#include <chrono>


using Complex = std::complex<double>;

std::vector<double> linspace(double start, double end, int steps) {
    std::vector<double> result;
    if (steps < 2) {
        // Edge case: if only one step, just return start
        result.push_back(start);
        return result;
    }
    
    result.reserve(steps);
    double step_size = (end - start) / (steps - 1);
    
    for (int i = 0; i < steps; ++i) {
        result.push_back(start + i * step_size);
    }
    
    return result;
}


void benchmark_biharm()
{
    double Lx = 1.10, Ly = 0.8, Lz = 5e-3;
    double E    = 9.0e+9 ;        // Young's mod [Pa]
    double nu   = 0.3 ;             // poisson's ratio
    double BC   = 1e15;   // elastic constants around the edges
    
    std::string filename = "/Users/admin/Documents/GitHub/magpie-cpp/data/benchmark_biharm_cpp.csv";
    std::vector<std::string> header = {"GridPoints", "ComputeTime"};
    
    std::ifstream infile(filename);
    bool file_exists = infile.good();
    infile.close();
    
    std::ofstream outfile(filename, std::ios::out);
    
    if (!file_exists) {
        outfile << header[0] << ',' << header[1] << "\n";
    }
    
    outfile.close();
    int Nxy = 0;
    
    //    for (auto val : linspace(0.005,0.02,1)) {
    for (auto val : linspace(0.05,0.001,2000)){
        double h    = std::sqrt(Lx*Ly)*val; // Grid Spacing
        
        int Nx = int(Lx / h);
        int Ny = int(Ly / h);
        
        if (Nx*Ny == Nxy)
            continue;
                        
        Nxy = Nx*Ny;
        
        double seconds = 1e15;
        
        for (int i = 0; i < 20; i++) {
            
            auto start = std::chrono::steady_clock::now();
            
            
            Eigen::SparseMatrix<double> biharm(Nxy,Nxy);
            bhmat(biharm,{BC,BC,BC,BC,BC,BC,BC,BC}, {Nx,Ny}, h, Lz, E, nu);
            
            
            auto end = std::chrono::steady_clock::now();
            double new_record  = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1.0e-9;
            std::cout << biharm.coeff(1,1);
                                    
            if (new_record < seconds) {
                seconds = new_record;
            }
        }
        
        std::cout<<'\n'<<'\n';
        std::cout << Nxy << ','  << ',' << seconds << '\n';
        outfile.open(filename, std::ios::app);
        outfile << (Nxy) << ','<< std::fixed << std::setprecision(15)<< seconds << '\n';
        outfile.close();
    }
}


void benchmark_eigs()
{
    double Lx = 1.10, Ly = 0.8, Lz = 5e-3;
    double E    = 9.0e+9 ;        // Young's mod [Pa]
    double nu   = 0.3 ;             // poisson's ratio
    double BC   = 1e15;   // elastic constants around the edges
    double h    = std::sqrt(Lx*Ly)*0.01; // Grid Spacing
    int Nx = int(Lx / h);
    int Ny = int(Ly / h);
    int Nxy = Nx*Ny;
  
    
    std::string filename = "/Users/admin/Documents/GitHub/magpie-cpp/data/benchmark_eigs_cpp_gen.csv";
    std::vector<std::string> header = {"Eigenvalues", "ComputeTime"};
    
    std::ifstream infile(filename);
    bool file_exists = infile.good();
    infile.close();
    
    std::ofstream outfile(filename, std::ios::out);
    
    if (!file_exists) {
        outfile << header[0] << ',' << header[1] << "\n";
    }
    
    outfile.close();
    
    Eigen::SparseMatrix<double> biharm(Nxy,Nxy);
    bhmat(biharm,{BC,BC,BC,BC,BC,BC,BC,BC}, {Nx,Ny}, h, Lz, E, nu);
    Spectra::SparseGenRealShiftSolve<double> op(biharm);
//    Spectra::SparseSymShiftSolve<double> op(biharm);
    
    
        for (int Nm = 10; Nm < 510; Nm+=10){
//    for (int Nm = 10; Nm < 20; Nm+=10){
        
        double seconds = 1e15;
        
        for (int i = 0; i < 20; i++){
//            auto start = std::chrono::steady_clock::now();
            auto start = std::chrono::high_resolution_clock::now();
            Spectra::GenEigsRealShiftSolver<Spectra::SparseGenRealShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0);
//            Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0);
            eigs.init();
            auto nconv = eigs.compute(Spectra::SortRule::LargestMagn);
            
//            auto end = std::chrono::steady_clock::now();
            auto end = std::chrono::high_resolution_clock::now();
            double new_record  = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1.0e-9;
                        
            if(eigs.info() == Spectra::CompInfo::Successful) // Retrieve results
            {
                auto evalues = eigs.eigenvalues();
                auto evectors = eigs.eigenvectors();
                std::cout <<nconv<<' ' << evalues[0];
                if (new_record < seconds) {
                    seconds = new_record;
                }
            }
        }
        
        outfile.open(filename, std::ios::app);
        std::cout << '\n' << '\n';
        std::cout << Nm << ',' << seconds << '\n';
        outfile << Nm << ','<< std::fixed << std::setprecision(15)<< seconds << '\n';
        outfile.close();
    }
}


double eig2hz (const double eig,double E, double rho, double Lz, double nu)
{
    double D = E * Lz*Lz*Lz / 12.0 / (1.0-(nu*nu));
    return std::sqrt(std::abs(eig))*std::sqrt(D/rho/Lz);
}

void simulate()
{
    double Lx = 1.10, Ly = 0.8, Lz = 5e-3;
    double E    = 9.0e+9 ;        // Young's mod [Pa]
    double nu   = 0.3 ;             // poisson's ratio
    double BC   = 1e15;   // elastic constants around the edges
    double rho  = 8765;
    double h    = std::sqrt(Lx*Ly)*0.01; // Grid Spacing
    int Nx = int(Lx / h);
    int Ny = int(Ly / h);
    int Nxy = Nx*Ny;
    int Nm = 500;
    
    Eigen::SparseMatrix<double> biharm(Nxy,Nxy);
    bhmat(biharm,{BC,0.0,0.0,0.0,BC,0.0,0.0,0.0}, {Nx,Ny}, h, Lz, E, nu);
    
    {
        std::string filename = "/Users/admin/Documents/GitHub/magpie-cpp/data/eigs_cpp_gen_cantilever.csv";
        
        Spectra::SparseGenRealShiftSolve<double> op(biharm);
        Spectra::GenEigsRealShiftSolver<Spectra::SparseGenRealShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0);
        
        eigs.init();
        auto nconv = eigs.compute(Spectra::SortRule::LargestMagn);
        if(eigs.info() == Spectra::CompInfo::Successful) // Retrieve results
        {
            auto evalues = eigs.eigenvalues();
            auto evectors = eigs.eigenvectors();
            
            std::ofstream outfile(filename, std::ios::out);
                                    
            for (int i = 0; i < nconv; i++) {
                outfile << std::fixed << std::setprecision(15) 
                << eig2hz(evalues[i].real(), E, rho, Lz, nu) << '\n';
            }
            
            
            outfile.close();
        }
    }
    
    
    {
        std::string filename = "/Users/admin/Documents/GitHub/magpie-cpp/data/eigs_cpp_sym_cantilever.csv";
        
        Spectra::SparseSymShiftSolve<double> op(biharm);
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0);

        
        eigs.init();
        auto nconv = eigs.compute(Spectra::SortRule::LargestMagn);
        if(eigs.info() == Spectra::CompInfo::Successful) // Retrieve results
        {
            auto evalues = eigs.eigenvalues();
            auto evectors = eigs.eigenvectors();
            
            std::ofstream outfile(filename, std::ios::out);
                                    
            for (int i = 0; i < nconv; i++) {
                outfile << std::fixed << std::setprecision(15)
                << eig2hz(evalues[i], E, rho, Lz, nu) << '\n';
            }
            
            
            outfile.close();
        }
    }
}





int main(int argc, const char * argv[])
{
    simulate();
//    benchmark_eigs();
//    benchmark_biharm();
    
//    double Lx = 1.10, Ly = 0.8, Lz = 5e-3;
//    double E    = 9.0e+9 ;        // Young's mod [Pa]
//    double nu   = 0.3 ;             // poisson's ratio
//    int Nm   = 400;               // number of modes to compute
//    double h    = std::sqrt(Lx*Ly)*0.01; // Grid Spacing
//    double BC   = 1e15;   // elastic constants around the edges
//    
//    int Nx = 110;
//    int Ny = 88;
//    
//    
//    Eigen::SparseMatrix<double> biharm(Ny*Nx,Ny*Nx);
//    bhmat(biharm,{BC,BC,BC,BC,BC,BC,BC,BC}, {Nx,Ny}, h, Lz, E, nu);
//    
//    Spectra::SparseGenRealShiftSolve<double> op(biharm);
//    Spectra::GenEigsRealShiftSolver<Spectra::SparseGenRealShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0);
    
    //    auto start = std::chrono::steady_clock::now();
//                    Spectra::SparseSymShiftSolve<double> op(biharm);
//                    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, Nm, (2*Nm) + 1,0.0);
    //
    //
    
    //
    //        //        switch (eigs.info()) {
    //        //            case Spectra::CompInfo::Successful:
    //        //                std::cout << "Success" << std::endl;
    //        //                break;
    //        //            case Spectra::CompInfo::NumericalIssue:
    //        //                std::cout << "NumericalIssue" << std::endl;
    //        //                break;
    //        //            case Spectra::CompInfo::NotConverging:
    //        //                std::cout << "NotConverging" << std::endl;
    //        //                break;
    //        //            case Spectra::CompInfo::NotComputed:
    //        //                std::cout << "NotComputed" << std::endl;
    //        //                break;
    //        //        }
    //
    
    ////    }
    //
    //    auto end = std::chrono::steady_clock::now();
    //    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //
    //    std::cout << "Elapsed time: " << elapsed_ms.count() << " ms\n";
    
    return 0;
}
