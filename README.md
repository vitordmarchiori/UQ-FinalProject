# Uncertainty Quantification in a Simple Damped Pendulum

This repository contains the final project for the course *Quantification of Uncertainties in Physical Systems*, June 2025.

<p align="center">
  <img src="Animations/pendulo_animacao.gif" alt="Oscillating Pendulum" width="100%">
</p>

The goal is to quantify the uncertainty propagation from an uncertain pendulum length to the dynamic response of the system.

---

## 📌 Description

The physical system consists of a damped simple pendulum with uncertain length \( L \sim \text{Gamma}(\nu, \theta) \). The governing equation is derived and solved for multiple realizations using the Monte Carlo method.

## 📊 Requirements

* MATLAB R2021a or newer
* Toolboxes: `Statistics and Machine Learning`, `ODE`, `MATLAB Base`

## 🔬 Topics Covered

* Linearized dynamics
* Damped oscillator
* Gamma-distributed uncertainty
* Monte Carlo simulation
* Propagation of uncertainty
* Kernel density estimation
* Visualization of confidence bands

---

## 📚 Acknowledgements

This project was developed under the supervision of Prof. Dr. Americo Barbosa da Cunha Junior, during the short course "Quantification of Uncertainties in Physical Systems" (2025).
