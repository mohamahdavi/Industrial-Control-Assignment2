# Industrial Control – Project Assignment 2

This repository contains my work for **Project Assignment 2** of the *Industrial Control* course (25791), instructed by **Behzad Ahi**, Electrical Engineering Department, Sharif University of Technology.  
The assignment covers a wide range of classical and modern control topics, with a strong emphasis on PID controller tuning, system identification, and practical implementation in MATLAB/Simulink.

## Contents
- **`HW2_Industrial_Control_Solutions.pdf`** – Full report with detailed derivations, explanations, and results (in Persian). All MATLAB code is appended at the end of the PDF.
- **MATLAB source codes (`*.m`)** – Separate script files for each question (or section), allowing direct reproduction of simulations and analyses.
- Additional datasets (e.g., `Ident.mat`, `GrowthData*.csv`) if required by the problems.

## Problems Covered

| # | Topic | Key Techniques |
|---|-------|----------------|
| 1 | PID Tuning using CHR Method | Step response, tracking, no overshoot vs. 20% overshoot |
| 2 | Simulating System Response & Model Fitting | Step response of high‑order system; first‑, second‑, and fourth‑order approximations |
| 3 | Method of Moments for System Identification | Swept‑frequency input, moment calculations, reduced‑order modelling |
| 4 | Reaction Curve Based PID Tuning | Ziegler‑Nichols & Cohen‑Coon, sensitivity analysis, stability margins, disturbance rejection |
| 5 | PID Synthesis via Transfer Function Matching | Converting general second‑order controller to standard PID form |
| 6 | Setpoint Tracking via Pole Placement | Augmented system, Sylvester matrix, ramp disturbance rejection |
| 7 | Windup and Anti‑Windup | Back‑calculation, conditional integration, Simulink PID block |
| 8 | MATLAB Identification Toolbox | Using `ident` to estimate transfer functions, goodness‑of‑fit |
| 9 | Model Order Reduction | Model Reducer app, frequency response comparison |
| 10 | Curve Fitting Toolbox | Biological growth models (Logistic, Gompertz, Richards, Bertalanffy) |
| 11 | PID Controller Optimization (GA) | Genetic Algorithm, cost function design, robustness analysis |
| 12 | Advanced PID Structures | Derivative on output, set‑point weighting, multi‑objective Pareto optimization |

## How to Use
1. Open the **PDF** for the full solution walk‑through (Persian).
2. Run the **MATLAB scripts** in the same directory to reproduce the figures and analyses.  
   *Make sure required toolboxes (Control System, Optimization, System Identification, Curve Fitting, etc.) are installed.*

## Course Details
- **Course:** Industrial Control (25791)  
- **Instructor:** Prof. Behzad Ahi  
- **Date:** May 2026 

## Author
- Mohammad Reza Mahdavi

*Status: Completed*
