#pragma once
#include "fem/Model.hpp"
#include <Eigen/Sparse>
#include <vector>

struct SolverResult {
    bool   converged = false;
    std::vector<double> displacements;  // displacement vector u [totalDOF]
    std::string message;
};

class Solver {
public:
    SolverResult solve(
        Eigen::SparseMatrix<double>& K,
        Eigen::VectorXd&             f,
        const std::vector<int>&      fixedDOFs
    );

private:
    void applyPenalization(
        Eigen::SparseMatrix<double>& K,
        Eigen::VectorXd&             f,
        const std::vector<int>&      fixedDOFs
    );
};