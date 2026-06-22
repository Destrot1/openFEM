#include "fem/Solver.hpp"
#include <Eigen/SparseLU>
#include <fmt/core.h>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <atomic>

void Solver::applyPenalization(
    Eigen::SparseMatrix<double>& K,
    Eigen::VectorXd&             f,
    const std::vector<int>&      fixedDOFs)
{
    // Find the maximum value on the diagonal
    double maxK = 0.0;
    for (int i = 0; i < K.rows(); ++i)
        maxK = std::max(maxK, std::abs(K.coeff(i, i)));

    double P = maxK * 1e12;
    if (P < 1.0) P = 1e12;  // fallback if K is all zeros

    fmt::print("[Solver] Penalization factor: {:.3e}\n", P);

    // Apply penalization -- K must be uncompressed to modify entries
    K.uncompress();
    for (int dof : fixedDOFs) {
        if (dof < 0 || dof >= K.rows()) continue;
        K.coeffRef(dof, dof) += P;
        f(dof) = 0.0;  // zero force on the constrained DOF
    }

    // TODO: penalization only supports homogeneous BCs (fixed = 0).
    // To impose a prescribed nonzero displacement d on a DOF, set
    // f(dof) = P * d instead of 0. Extend the signature accordingly.
}

SolverResult Solver::solve(
    Eigen::SparseMatrix<double>& K,
    Eigen::VectorXd&             f,
    const std::vector<int>&      fixedDOFs)
{
    SolverResult result;

    if (K.rows() != K.cols()) {
        result.message = "K is not square";
        return result;
    }
    if (K.rows() != f.size()) {
        result.message = "K and f size mismatch";
        return result;
    }

    fmt::print("[Solver] System size: {}x{}\n", K.rows(), K.cols());
    fmt::print("[Solver] Fixed DOFs: {}\n", fixedDOFs.size());

    // -- Apply boundary conditions ---------------------------------
    applyPenalization(K, f, fixedDOFs);

    // -- Solve K*u = f with SparseLU -------------------------------
    // SparseLU is a direct solver -- exact, not iterative.
    // Suitable for systems up to ~100k DOF.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    auto t0 = std::chrono::steady_clock::now();
    auto elapsedSec = [](auto since) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
    };

    fmt::print("[Solver] Analyzing sparsity pattern...\n");
    solver.analyzePattern(K);
    fmt::print("[Solver] Pattern analyzed ({:.1f}s)\n", elapsedSec(t0));

    // Factorization is the slow, silent part for a large system, with
    // no way to get incremental progress out of Eigen's SparseLU --
    // it's one blocking call. Run it on a separate thread so the main
    // thread can print a heartbeat in the meantime, just so a long
    // wait doesn't look like the program has frozen.
    auto t1 = std::chrono::steady_clock::now();
    fmt::print("[Solver] Factorizing {}x{} ({} non-zeros) -- this can take a "
               "while for a large system...\n", K.rows(), K.cols(), K.nonZeros());

    std::atomic<bool> factorizeDone{false};
    std::thread heartbeat([&]() {
        while (!factorizeDone.load()) {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            if (!factorizeDone.load())
                fmt::print("[Solver] ... still factorizing ({:.0f}s elapsed)\n", elapsedSec(t1));
        }
    });
    solver.factorize(K);
    factorizeDone = true;
    heartbeat.join();

    if (solver.info() != Eigen::Success) {
        result.message = "Factorization failed -- matrix may be singular";
        fmt::print("[Solver] ERROR: {}\n", result.message);
        return result;
    }
    fmt::print("[Solver] Factorized ({:.1f}s)\n", elapsedSec(t1));

    auto t2 = std::chrono::steady_clock::now();
    fmt::print("[Solver] Back-substituting...\n");
    Eigen::VectorXd u = solver.solve(f);

    if (solver.info() != Eigen::Success) {
        result.message = "Solve failed";
        fmt::print("[Solver] ERROR: {}\n", result.message);
        return result;
    }
    fmt::print("[Solver] Done ({:.1f}s)\n", elapsedSec(t2));

    // -- Copy results ----------------------------------------------
    result.converged    = true;
    result.displacements.assign(u.data(), u.data() + u.size());
    result.message      = "OK";

    fmt::print("[Solver] Solution converged.\n");
    fmt::print("[Solver] Max displacement: {:.6e} mm\n", u.cwiseAbs().maxCoeff());

    return result;
}