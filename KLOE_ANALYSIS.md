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
> **BDT test mode:** `root -l -q $HOME/run_bdt/test_bdt.C #Outputs: $HOME/Desktop/input_bdt_TDATA_norm/`
> **BDT analysis:** `root - l -q $HOME/input_bdt.sh #Output: $HOME/Destop/KLOE_BDT/plots_bdt/` 


## 📐 Analysis Flow

### 1.Data Preparation
#### Input Raw Data 
```bash
script/listpath_chain.sh # listing path of raw data root files stored as a text input file. Or listpath_norm.sh for larger data samples
# Outputs: $HOME/Destop/KLOE_BDT/path_norm/*
```

#### Create ROOT Files 
```bash
root -l -b -q run_bdt/Process.C #prompt, small samples
# Outputs: $HOME/Desktop/sig.root

./run_bdt/script/input_bdt.sh 
# Outputs: $HOME/Desktop/input_bdt_TDATA_chain/*
```

