#include <ControllerToolbox.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

class BoilerTurbine
{
public:
    const float Ts = 1.0f; // Sampling time in seconds

    // Commanded valve positions (set externally each step)
    float u1 = 0.5f, u2 = 0.5f, u3 = 0.5f;

    // Plant states and outputs
    float x1, x2, x3;
    float y1, y2, y3;

    // Rate-of-change of valve positions (read-back after constrain_valve_rate)
    float du1 = 0.0f, du2 = 0.0f, du3 = 0.0f;

    inline void constrain_valve()
    {
        u1 = std::clamp(u1, 0.0f, 1.0f);
        u2 = std::clamp(u2, 0.0f, 1.0f);
        u3 = std::clamp(u3, 0.0f, 1.0f);
    }

    // Computes du = u - u_prev, clamps the rate, then applies the
    // rate-limited increment so u reflects the actual achievable position.
    // du1/du2/du3 are updated as diagnostic outputs.
    inline void constrain_valve_rate()
    {
        du1 = std::clamp(u1 - u1_prev_, -0.007f, 0.007f);
        du2 = std::clamp(u2 - u2_prev_, -2.0f, 0.02f);
        du3 = std::clamp(u3 - u3_prev_, -0.05f, 0.05f);
        u1 = u1_prev_ + du1;
        u2 = u2_prev_ + du2;
        u3 = u3_prev_ + du3;
        u1_prev_ = u1;
        u2_prev_ = u2;
        u3_prev_ = u3;
    }

    void update()
    {
        float x1_98 = std::pow(x1, 9.0f / 8.0f);

        float dx1 = -0.0018f * u2 * x1_98 + 0.9f * u1 - 0.15f * u3;
        float dx2 = (0.073f * u2 - 0.016f) * x1_98 - 0.1f * x2;
        float dx3 = (141.0f * u3 - (1.1f * u2 - 0.19f) * x1) / 85.0f;

        // Outputs from current state (all at time k, before the Euler step)
        y1 = x1;
        y2 = x2;
        float acs = ((1.0f - 0.001538f * x3) * 0.8f * x1 - 25.6f) / (x3 * (1.0394f - 0.0012304f * x1));
        float qe = (0.854f * u2 - 0.147f) * x1 + 45.59f * u1 - 2.514f * u3 - 2.096f;
        y3 = 0.05f * (0.13073f * x3 + 100.0f * acs + qe / 9.0f - 67.975f);

        // Forward-Euler state update: x[k+1] = x[k] + Ts * dx
        x1 += Ts * dx1;
        x2 += Ts * dx2;
        x3 += Ts * dx3;
    }

private:
    float u1_prev_ = 0.5f, u2_prev_ = 0.5f, u3_prev_ = 0.5f;
};

void random_valve(float &u1, float &u2, float &u3)
{
    // Random between change rate
    float du1 = ((float)rand() / RAND_MAX) * 0.014f - 0.007f; // [-0.007, 0.007]
    float du2 = ((float)rand() / RAND_MAX) * 0.04f - 0.02f;   // [-0.02, 0.02]
    float du3 = ((float)rand() / RAND_MAX) * 0.1f - 0.05f;    // [-0.05, 0.05]

    u1 += du1;
    u2 += du2;
    u3 += du3;
}

void plant_model_data()
{
    BoilerTurbine bt;

    // Initial conditions
    bt.x1 = 100.0f; // Drum pressure
    bt.x2 = 50.0f;  // Electric power
    bt.x3 = 20.0f;  // Water level

    // Control inputs
    bt.u1 = 0.5f; // Fuel flow valve position
    bt.u2 = 0.5f; // Steam control valve position
    bt.u3 = 0.5f; // Feedwater flow valve position

    std::ofstream data_file("boiler_turbine_data.csv");
    data_file << "Time,Drum Pressure,Electric Power,Water Level Deviation,u1,u2,u3\n";

    for (int i = 0; i < 100; ++i)
    {
        random_valve(bt.u1, bt.u2, bt.u3);
        bt.constrain_valve_rate();
        bt.constrain_valve();
        bt.update();
        data_file << i * bt.Ts << "," << bt.y1 << "," << bt.y2 << "," << bt.y3 << "," << bt.u1 << "," << bt.u2 << "," << bt.u3 << "\n";
    }
}

int main()
{
    plant_model_data();
    std::cout << "Boiler-turbine data generated in boiler_turbine_data.csv\n";
    return 0;
}