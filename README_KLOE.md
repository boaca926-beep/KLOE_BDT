## 📚 Table of Contents
- [Overview](#-overview)
- [Standard Analysis ($\chi^{2}$)](#-standard-analysis-χ)
- [BDT Analysis](#-bdt-analysis)
- [Data Preparation](#-data-preparation)
- [BDT Training & Evaluation](#-bdt-training--evaluation)
- [Performance Metrics](#-performance-metrics)
- [Quick Start](#-quick-start)
- [Repository Structure](#-repository-structure)

## 💡 Overview
This project analyzes the $e^{+}e^{-}\to\pi^{+}\pi^{-}\pi^{0}\gamma$ ISR process, using a $1.7~\text{fb}^{-1}$ data sample collected at KLOE.

<div align="center">
<img src="plots_ref/data_flow_new_IV.png" width="500" alt="ROC & AUC"/>
<br/>
<em>Figure 1: Analysis workflow for Monte Carlo (MC) simulated events and data, beginning with events passing the trigger and selected by the FILFO and KSL streams. The final state of all events is reconstructed using a kinematic fit that includes e⁺e⁻ → π⁺π⁻π⁰ photon pairing.</em>

</div>

## 🚀 Quick Start {#-quick-start}
> **Quick start:** `uv sync && uv run main_initialize_kloe_opti.py`

## 📐 Analysis Flow


