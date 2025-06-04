# 🎯 Uncertainty Quantification in a Simple Damped Pendulum

This repository contains the final project for the course *Quantification of Uncertainties in Physical Systems*, conducted in June 2025.

<p align="center">
  <img src="Animations/pendulo_animacao.gif" alt="Oscillating Pendulum" width="100%">
</p>

The objective is to model a damped simple pendulum and quantify the effects of uncertainty in its length through a probabilistic framework and Monte Carlo simulation.

---

## 📌 Description

A simple damped pendulum is analyzed, considering uncertainty in its length \( L \), modeled as a Gamma-distributed random variable. The angular response \( \theta(t) \) is obtained via numerical integration of the linearized equation of motion, and statistical quantities of interest are computed from the ensemble of trajectories.

---

## 📊 Topics Covered

* Analytical modeling of a damped pendulum
* Small-angle approximation and linearization
* Uncertainty modeling using Gamma distributions
* Monte Carlo simulation (MC)
* Propagation of uncertainty in dynamic systems
* Statistical post-processing (mean, std, confidence intervals)
* Histogram and kernel density estimation (KDE)

---

## 🧰 Requirements

* MATLAB R2021a or newer
* Toolboxes:

  * Statistics and Machine Learning
  * ODE Suite
  * MATLAB Base

---

## 📚 Acknowledgements

This project was developed as the final assignment for the short course:

> **Quantification of Uncertainties in Physical Systems**
> Instructor: Prof. Dr. Americo Barbosa da Cunha Junior
> Institution: São Paulo State University (UNESP), School of Engineering, Ilha Solteira

---

## 🔗 Link to Report

The full technical report (in LaTeX) is available in the `/Code and Report` folder and will be submitted alongside this repository as part of the course evaluation.
