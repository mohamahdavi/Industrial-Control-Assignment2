# PID Tuning, System Identification & Multi-Objective Optimization

[![MATLAB](https://img.shields.io/badge/MATLAB-Simulink-0076a8.svg)](#getting-started)
[![Report PDF](https://img.shields.io/badge/Report-Technical%20PDF-b31b1b.svg)](HW2_401102637-Mahdavi/HW2_Industrial_Control.pdf)
[![Course](https://img.shields.io/badge/Course-Industrial%20Control%20(SUT)-1f4e79.svg)](#author)

End-to-end industrial control design and system identification suite covering empirical PID tuning rules, experimental model identification, model-order reduction, Simulink anti-windup implementation, and multi-objective Genetic Algorithm (GA) Pareto optimization. Developed for the **Industrial Control** course at **Sharif University of Technology**.

> [!IMPORTANT]
> **Design Objective:** Bridge classical industrial tuning heuristics with modern identification and multi-objective evolutionary optimization while addressing practical actuator saturation via Simulink anti-windup architectures.

---

## Key Components & Technical Highlights

| Module | Methods & Artifacts | Key Details |
| :--- | :--- | :--- |
| **Classical PID Tuning** | CHR, Ziegler–Nichols & Cohen–Coon | Comparative step-response benchmarking across empirical industrial tuning rules (`Q1.m`–`Q12.m`). |
| **System Identification** | Method of Moments & Identification Toolbox | Experimental parameter estimation from time-series datasets (`GrowthData1.csv`–`GrowthData4.csv`, `Ident.mat`) saved in `estimated_model.mat`, followed by model-order reduction. |
| **Industrial PID Structures** | Setpoint Weighting & Derivative-on-Output | Pole-placement setpoint tracking paired with structural modifications to prevent derivative kick and overshoot. |
| **Actuator Saturation & Anti-Windup** | Simulink Implementation | Closed-loop validation of back-calculation and conditional integration anti-windup schemes in [`PID_antiwindup_simulink.slx`](HW2_401102637-Mahdavi/SimulationFiles/PID_antiwindup_simulink.slx). |
| **Multi-Objective Optimization** | Genetic Algorithm (Pareto Front) | Trade-off optimization of PID gains balancing tracking speed against control effort (`pareto_results.mat`). |

---

## Scope & Limitations

- **Simulation Environment:** Evaluated in continuous/discrete-time MATLAB scripts (`Q1.m`–`Q12.m`) and Simulink (`PID_antiwindup_simulink.slx`) using provided experimental datasets (`GrowthData1`–`4`).
- **Baseline Benchmarking:** Classical tuning heuristics serve as baselines to demonstrate the disturbance-rejection and saturation-handling improvements of anti-windup and Pareto-optimal designs.

---

## Getting Started

```bash
git clone https://github.com/mohamahdavi/Industrial-Control-Assignment2.git
cd Industrial-Control-Assignment2/HW2_401102637-Mahdavi/SimulationFiles

# Open MATLAB and run the problem scripts Q1.m through Q12.m
# Open PID_antiwindup_simulink.slx in Simulink for anti-windup simulations
```

---

## Repository Structure

```text
.
├── README.md
└── HW2_401102637-Mahdavi/
    ├── HW2_Industrial_Control.pdf        # Compiled technical report
    └── SimulationFiles/
        ├── Q1.m ... Q12.m                # MATLAB scripts for tuning, ID, reduction & GA optimization
        ├── PID_antiwindup_simulink.slx   # Simulink model for actuator saturation & anti-windup
        ├── GrowthData1.csv ... 4.csv     # Experimental step/growth datasets for system identification
        ├── Ident.mat                     # Identification workspace data
        ├── estimated_model.mat           # Identified & reduced-order transfer function models
        └── pareto_results.mat            # Multi-objective GA Pareto-optimal PID solutions
```

---

## Author

- **Mohammadreza Mahdavi** — B.Sc. Student in Electrical Engineering (Control Systems), Sharif University of Technology
