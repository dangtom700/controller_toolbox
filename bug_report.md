# Controller Toolbox Improvement Plan & Bug Report

This document outlines the critical gaps currently existing in the C++ Controller Toolbox. The gaps are divided into two categories: **Control Theory & Analysis Tools** (missing methods to assess stability and performance) and **Software Architecture Limitations** (problems preventing the code from being deployed in hard real-time environments).

---

## Part 1: Missing Control Theory & Analysis Tools

While the toolbox provides excellent controller implementations (PID, MPC, LQG, ADRC) and time-domain simulation capabilities, it severely lacks the analytical tools required to *prove* stability and response speed without running brute-force simulations.

### 1. Stability Margin Calculators (Frequency Domain)
* **The Problem:** There is currently no way to calculate the Gain Margin (GM) and Phase Margin (PM) of a system. 
* **The Fix:** Implement Bode plot and Nyquist stability criterion analyzers. Add methods to `StateSpace` and `TransferFunction` classes that automatically compute the crossover frequencies and return the stability margins.

### 2. Automated Time-Domain Metric Extraction
* **The Problem:** In our simulation scripts, we had to manually write `for`-loops to track the maximum value and calculate the Integral Squared Error (ISE).
* **The Fix:** Implement a `MetricsAnalyzer` class that takes raw step-response data (arrays of time and output) and automatically calculates:
  * Rise Time ($t_r$)
  * Settling Time ($t_s$)
  * Peak Overshoot ($M_p$)
  * Steady-State Error ($e_{ss}$)

### 3. Pole/Zero Extraction and Root Locus Analysis
* **The Problem:** The toolbox can construct `StateSpace` models, but it does not expose an easy way to extract the eigenvalues (poles) of the `A` matrix to check for stability mathematically.
* **The Fix:** Add a `getPoles()` and `getZeros()` method to all plant models. Implement a function that verifies if all discrete poles lie strictly inside the unit circle ($|z| < 1$).

### 4. Robust Sensitivity Analysis ($H_\infty$ norms)
* **The Problem:** Evaluating robustness currently requires running Monte Carlo simulations (as done in `ex22_full_pipeline_robustness.cpp`).
* **The Fix:** Implement solvers to calculate the peak of the Sensitivity ($S$) and Complementary Sensitivity ($T$) functions (i.e., calculating the $H_\infty$ norm of the closed-loop system) to provide a strict mathematical guarantee of robustness against modeling errors.

### 5. Discrete Lyapunov Equation Solver
* **The Problem:** No strict mathematical proof of stability for MIMO LQR/LQG systems in the presence of perturbations.
* **The Fix:** Provide an algebraic solver for the discrete Lyapunov equation ($A^T P A - P + Q = 0$).

---

## Part 2: Software Architecture Problems & Real-Time Constraints

Even if the mathematical models are perfect, the current C++ implementation cannot safely be deployed onto an embedded microcontroller or a hard real-time operating system (RTOS) due to several architectural software flaws.

### 1. Non-Deterministic Memory Allocation
* **The Problem:** The toolbox relies heavily on dynamic memory allocation during runtime. `Eigen::MatrixXd` and `std::vector` resize operations under the hood (e.g., inside the MPC or LQR solvers) call `malloc` or `new`. In a real-time loop, memory allocation can cause unpredictable latency spikes, leading to missed control deadlines.
* **The Fix:** Enforce "Zero-Allocation" execution. The toolbox must pre-allocate all matrices and buffers at startup using fixed-size types (e.g., `Eigen::Matrix<double, N, M>`) or strictly prevent resizing inside the `compute()` methods.

### 2. Lack of Thread-Safety and Lock-Free Updates
* **The Problem:** Updating controller parameters on-the-fly (like the gain scheduling or adaptive MPC updates) is not thread-safe. If a background thread updates the `MPCParams` while the real-time thread is executing `compute()`, a race condition will corrupt the matrices.
* **The Fix:** Implement lock-free atomic pointers or a double-buffering mechanism so that background telemetry can inject new parameters without acquiring blocking mutexes that stall the control loop.

### 3. Missing Hardware Abstraction Layer (HAL)
* **The Problem:** The toolbox is purely mathematical. There is no standard interface to connect the controllers to physical ADCs (Analog-to-Digital Converters) or DAC/PWM actuators.
* **The Fix:** Define standard `ISensor` and `IActuator` C++ interfaces. This allows the exact same control logic to run against a simulated mathematical plant or a real EtherCAT/CAN-bus hardware device without altering the core logic.

### 4. Poor Exception Handling and No "Safe-State" Fallbacks
* **The Problem:** If a sensor returns `NaN` or a matrix inversion fails inside `DiscreteMPC`, the code currently does nothing to protect the hardware. It will likely pass the `NaN` to the actuator or throw an unhandled C++ exception, crashing the entire process.
* **The Fix:** Add rigorous `std::isnan` checks at all inputs. Implement a deterministic fallback state (e.g., freezing the actuator or switching to a highly conservative PID controller) when the advanced solver fails or detects anomalous data.

### 5. Blocking Telemetry and Logging
* **The Problem:** Writing step-response data to CSV files (like in `simulate_all.cpp`) uses `std::ofstream`, which relies on blocking I/O calls to the hard drive. Doing this inside a control loop destroys real-time performance.
* **The Fix:** Build a lock-free ring-buffer logging architecture. The real-time thread pushes data into memory in $O(1)$ time, and a separate low-priority thread pulls from the buffer to write to the disk or stream over a network.
