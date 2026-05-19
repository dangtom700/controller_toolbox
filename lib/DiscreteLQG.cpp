#include "DiscreteLQG.h"

namespace ctrl {

DiscreteLQG::DiscreteLQG(const StateSpace&       plant,
                         const LQRParams&        lqr_p,
                         const Eigen::MatrixXd&  Q_noise,
                         const Eigen::MatrixXd&  R_noise,
                         const Eigen::MatrixXd&  P0)
    : plant_(plant)
{
    lqr_ = std::make_unique<DiscreteLQR>(plant, lqr_p);
    kf_  = std::make_unique<KalmanFilter>(plant, Q_noise, R_noise, P0);
    x_ref_  = Eigen::VectorXd::Zero(plant.stateSize());
    u_prev_ = Eigen::VectorXd::Zero(plant.inputSize());
}

// Full step: KF predict → KF update → LQR control.
Eigen::VectorXd DiscreteLQG::step(const Eigen::VectorXd& y,
                                   const Eigen::VectorXd& u_prev,
                                   const Eigen::VectorXd& x_ref)
{
    // 1. Kalman predict using previous control
    kf_->predict(u_prev);

    // 2. Kalman update using current measurement
    kf_->update(y, u_prev);

    // 3. LQR feedback on estimated state
    const Eigen::VectorXd& xhat = kf_->state();
    Eigen::VectorXd ref = x_ref.size() == plant_.stateSize() ? x_ref : x_ref_;

    Eigen::VectorXd u = lqr_->compute(xhat, ref);

    // Store for internal convenience
    u_prev_ = u;
    return u;
}

// IController-compatible SISO wrapper.
// Requires setReference() and setUPrev() to be called first if needed.
double DiscreteLQG::compute(double y_scalar)
{
    Eigen::VectorXd y(plant_.outputSize());
    y.fill(y_scalar);
    Eigen::VectorXd u = step(y, u_prev_, x_ref_);
    return u(0);
}

void DiscreteLQG::reset()
{
    kf_->reset();
    u_prev_.setZero();
    x_ref_.setZero();
}

} // namespace ctrl
