#include "SystemAnalysis.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ctrl {

std::vector<std::complex<double>> SystemAnalysis::getPoles(const StateSpace& sys) {
    Eigen::EigenSolver<Eigen::MatrixXd> es(sys.A);
    Eigen::VectorXcd eigs = es.eigenvalues();
    std::vector<std::complex<double>> poles(eigs.size());
    for (int i = 0; i < eigs.size(); ++i) {
        poles[i] = eigs[i];
    }
    return poles;
}

bool SystemAnalysis::isDiscreteStable(const StateSpace& sys) {
    auto poles = getPoles(sys);
    for (const auto& p : poles) {
        if (std::abs(p) >= 1.0) {
            return false;
        }
    }
    return true;
}

Eigen::MatrixXd SystemAnalysis::solveDiscreteLyapunov(const Eigen::MatrixXd& A,
                                                      const Eigen::MatrixXd& Q) {
    // Solves A * P * A^T - P + Q = 0 by vectorizing:
    // (I - A \otimes A) vec(P) = vec(Q)
    int n = A.rows();
    if (n != A.cols() || n != Q.rows() || n != Q.cols()) {
        throw std::invalid_argument("solveDiscreteLyapunov: Matrix dimensions mismatch.");
    }
    
    // Kronecker product: A \otimes A
    Eigen::MatrixXd kronAA(n * n, n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            kronAA.block(i * n, j * n, n, n) = A(i, j) * A;
        }
    }
    
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n * n, n * n);
    Eigen::MatrixXd M = I - kronAA;
    
    // vec(Q)
    Eigen::VectorXd vecQ(n * n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            vecQ(j * n + i) = Q(i, j);
        }
    }
    
    // Solve M * vec(P) = vec(Q)
    Eigen::VectorXd vecP = M.colPivHouseholderQr().solve(vecQ);
    
    // Reshape back to matrix P
    Eigen::MatrixXd P(n, n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            P(i, j) = vecP(j * n + i);
        }
    }
    return P;
}

std::vector<std::complex<double>> SystemAnalysis::getFrequencyResponse(const StateSpace& sys,
                                                                       const std::vector<double>& freqs) {
    std::vector<std::complex<double>> response(freqs.size());
    int n = sys.A.rows();
    Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(n, n);
    
    // Convert to complex matrices
    Eigen::MatrixXcd Ac = sys.A.cast<std::complex<double>>();
    Eigen::MatrixXcd Bc = sys.B.cast<std::complex<double>>();
    Eigen::MatrixXcd Cc = sys.C.cast<std::complex<double>>();
    Eigen::MatrixXcd Dc = sys.D.cast<std::complex<double>>();
    
    for (size_t i = 0; i < freqs.size(); ++i) {
        // z = e^(j * w * Ts)
        std::complex<double> z = std::polar(1.0, freqs[i] * sys.Ts);
        
        // G(z) = C * (zI - A)^-1 * B + D
        Eigen::MatrixXcd resolvent = (z * I - Ac).inverse();
        Eigen::MatrixXcd G = Cc * resolvent * Bc + Dc;
        
        // For SISO
        if (G.rows() == 1 && G.cols() == 1) {
            response[i] = G(0, 0);
        } else {
            // For MIMO, just return the (0,0) element or handle differently
            response[i] = G(0, 0);
        }
    }
    return response;
}

StabilityMargins SystemAnalysis::calculateMargins(const StateSpace& sys) {
    StabilityMargins margins = {0.0, 0.0, 0.0, 0.0};
    
    // Logarithmic frequency grid from 10^-3 to Nyquist
    double w_min = 1e-3;
    double w_max = M_PI / sys.Ts;
    int num_points = 10000;
    
    std::vector<double> freqs(num_points);
    for(int i = 0; i < num_points; ++i) {
        freqs[i] = w_min * std::pow(w_max / w_min, static_cast<double>(i) / (num_points - 1));
    }
    
    auto resp = getFrequencyResponse(sys, freqs);
    
    double min_gain_dist = 1e9;
    double min_phase_dist = 1e9;
    
    for(size_t i = 0; i < resp.size(); ++i) {
        double mag = std::abs(resp[i]);
        double phase = std::arg(resp[i]) * 180.0 / M_PI; // degrees
        
        // Wrap phase to [-360, 0] assuming it starts near 0 and goes negative
        while(phase > 0) phase -= 360.0;
        
        // Find phase crossover: phase goes through -180
        if (std::abs(phase - (-180.0)) < min_phase_dist) {
            min_phase_dist = std::abs(phase - (-180.0));
            margins.wCrossoverGain = freqs[i];
            if (mag > 0) {
                margins.gainMarginDb = -20.0 * std::log10(mag);
            }
        }
        
        // Find gain crossover: mag goes through 1.0 (0 dB)
        if (std::abs(mag - 1.0) < min_gain_dist) {
            min_gain_dist = std::abs(mag - 1.0);
            margins.wCrossoverPhase = freqs[i];
            margins.phaseMarginDeg = 180.0 + phase;
        }
    }
    
    // If no crossover found, margins can be considered infinite
    if (min_phase_dist > 50.0) {
        margins.gainMarginDb = std::numeric_limits<double>::infinity();
    }
    if (min_gain_dist > 0.5) {
        margins.phaseMarginDeg = std::numeric_limits<double>::infinity();
    }
    
    return margins;
}

double SystemAnalysis::calculateHInfinityNorm(const StateSpace& sys) {
    // Grid-based search for peak singular value over frequency
    double w_min = 1e-3;
    double w_max = M_PI / sys.Ts;
    int num_points = 2000;
    
    int n = sys.A.rows();
    Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(n, n);
    Eigen::MatrixXcd Ac = sys.A.cast<std::complex<double>>();
    Eigen::MatrixXcd Bc = sys.B.cast<std::complex<double>>();
    Eigen::MatrixXcd Cc = sys.C.cast<std::complex<double>>();
    Eigen::MatrixXcd Dc = sys.D.cast<std::complex<double>>();
    
    double peak_mag = 0.0;
    
    for(int i = 0; i < num_points; ++i) {
        double w = w_min * std::pow(w_max / w_min, static_cast<double>(i) / (num_points - 1));
        std::complex<double> z = std::polar(1.0, w * sys.Ts);
        
        Eigen::MatrixXcd resolvent = (z * I - Ac).inverse();
        Eigen::MatrixXcd G = Cc * resolvent * Bc + Dc;
        
        // Max singular value is the induced 2-norm of the matrix
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(G);
        double max_sv = svd.singularValues()(0);
        
        if (max_sv > peak_mag) {
            peak_mag = max_sv;
        }
    }
    return peak_mag;
}

} // namespace ctrl
