#include "DiscretePID.h"
#include <algorithm>
#include <cmath>

namespace ctrl
{

    DiscretePID::DiscretePID(const PIDParams &params, double sampleTime)
        : p_(params), Ts_(sampleTime)
    {
        reset();
    }

    // Discrete PID — backward Euler integral + backward-Euler filtered derivative.
    //
    // Derivative filter recurrence (derived from C_D(s) = Kd·N·s/(s+N), backward Euler):
    //   α     = 1 / (1 + N·Ts)
    //   d[k]  = α·d[k-1] + Kd·N·α·(e[k] − e[k-1])
    //
    // Anti-windup back-calculation (Åström & Wittenmark):
    //   I[k+1] = I[k] + Ki·Ts·e[k] + Kb·(u_sat[k] − u_unsat[k])
    double DiscretePID::compute(double error)
    {
        // Hold last output on bad measurement rather than corrupting integral state
        if (!std::isfinite(error))
            return u_prev_;

        // Filtered derivative
        const double alpha = 1.0 / (1.0 + p_.N * Ts_);
        const double d_new = alpha * deriv_ + p_.Kd * p_.N * alpha * (error - e_prev_);

        // Unsaturated output (uses integral from previous step)
        const double u_unsat = p_.Kp * error + integral_ + d_new;

        // Output saturation
        const double u_sat = std::max(p_.uMin, std::min(p_.uMax, u_unsat));

        // Integral update with anti-windup back-calculation
        integral_ += p_.Ki * Ts_ * error + p_.Kb * (u_sat - u_unsat);

        // State updates
        deriv_ = d_new;
        e_prev_ = error;
        u_prev_ = u_sat;

        return u_sat;
    }

    void DiscretePID::reset()
    {
        integral_ = 0.0;
        deriv_ = 0.0;
        e_prev_ = 0.0;
        u_prev_ = 0.0;
    }

} // namespace ctrl
