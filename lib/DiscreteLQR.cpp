#include "DiscreteLQR.h"
#include <stdexcept>
#include <string>

namespace ctrl
{

    // ---------------------------------------------------------------------------
    // Value-iteration DARE solver
    //   P_{t+1} = A'P_t A − (A'P_t B)(R + B'P_t B)⁻¹(B'P_t A) + Q
    // Converges to the stabilising solution P∞ for stabilisable (A,B) + detectable (A,Q½).
    // Ref: Bertsekas "Dynamic Programming and Optimal Control" Vol 1, §1.3.
    // ---------------------------------------------------------------------------
    Eigen::MatrixXd DiscreteLQR::solveDARE(const Eigen::MatrixXd &A,
                                           const Eigen::MatrixXd &B,
                                           const Eigen::MatrixXd &Q,
                                           const Eigen::MatrixXd &R)
    {
        const int maxIter = 10000;
        const double tol = 1e-12;

        Eigen::MatrixXd P = Q;

        for (int iter = 0; iter < maxIter; ++iter)
        {
            const Eigen::MatrixXd AtP = A.transpose() * P;
            const Eigen::MatrixXd S = R + B.transpose() * P * B;
            const Eigen::MatrixXd Gain = S.ldlt().solve(B.transpose() * P * A);
            const Eigen::MatrixXd P_new = AtP * A - AtP * B * Gain + Q;

            const double err = (P_new - P).norm() / (1.0 + P.norm());
            P = P_new;

            if (err < tol)
                return P;
        }

        throw std::runtime_error(
            "DiscreteLQR: DARE value iteration did not converge in " +
            std::to_string(maxIter) + " iterations. "
                                      "Verify that (A,B) is stabilisable and (A,sqrt(Q)) is detectable.");
    }

    DiscreteLQR::DiscreteLQR(const StateSpace &plant, const LQRParams &params)
        : Ts_(plant.Ts), n_(plant.stateSize()), m_(plant.inputSize())
    {
        P_ = solveDARE(plant.A, plant.B, params.Q, params.R);
        const Eigen::MatrixXd S = params.R + plant.B.transpose() * P_ * plant.B;
        K_ = S.ldlt().solve(plant.B.transpose() * P_ * plant.A);
    }

    // u[k] = −K*(x − x_ref) + u_ff
    Eigen::VectorXd DiscreteLQR::compute(const Eigen::VectorXd &x,
                                         const Eigen::VectorXd &x_ref,
                                         const Eigen::VectorXd &u_ff) const
    {
        Eigen::VectorXd xe = x;
        if (x_ref.size() == n_)
            xe -= x_ref;

        Eigen::VectorXd u = -K_ * xe;
        if (u_ff.size() == m_)
            u += u_ff;
        return u;
    }

} // namespace ctrl
