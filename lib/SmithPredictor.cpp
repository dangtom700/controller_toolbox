#include "SmithPredictor.h"

namespace ctrl {

SmithPredictor::SmithPredictor(std::shared_ptr<IController> inner,
                               const StateSpace&             delayModel,
                               int                           delaySteps)
    : inner_(std::move(inner))
    , model_(delayModel)
    , d_(delaySteps)
    , Ts_(delayModel.Ts)
{
    x_model_ = Eigen::VectorXd::Zero(model_.stateSize());
    // Pre-fill the buffer with zeros (plant starts at rest)
    y_buf_.assign(d_, 0.0);
}

// Modified error:  e_sp = (r−y) + (ŷ_now − ŷ_delayed)
//                        = error + ŷ_now − ŷ_now[k−d]
//
// The correction term removes the delay from the error signal presented to
// the inner controller, giving it an effectively delay-free closed loop.
double SmithPredictor::compute(double error)
{
    // Current model output ŷ[k] from the delay-free plant model
    const double y_now     = (model_.C * x_model_ + model_.D *
                               Eigen::VectorXd::Zero(model_.inputSize()))(0);
    // d-step delayed model output (front of circular buffer)
    const double y_delayed = y_buf_.front();

    // Modified error: routes inner controller as if there were no dead-time
    const double e_sp = error + (y_now - y_delayed);

    // Compute control action via the inner controller
    const double u = inner_->compute(e_sp);

    // Advance internal model state: x̂[k+1] = A·x̂[k] + B·u[k]
    Eigen::VectorXd uv(model_.inputSize());
    uv.fill(u);
    x_model_ = model_.A * x_model_ + model_.B * uv;

    // Shift delay buffer: discard oldest, append current ŷ[k]
    y_buf_.pop_front();
    y_buf_.push_back(y_now);

    return u;
}

void SmithPredictor::reset()
{
    inner_->reset();
    x_model_.setZero();
    y_buf_.assign(d_, 0.0);
}

} // namespace ctrl
