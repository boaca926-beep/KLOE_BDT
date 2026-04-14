# 📄 KLOE Analysis (Standard)
## 💡 Description
- ROOT framework based analysis

- Study of $e^{+}e^{-}\to\pi^{+}\pi^{-}\pi^{0}\gamma\to\pi^{+}\pi^{-}\gamma\gamma\gamma$ ISR process.

### $\pi^{0}$ Reconstruction using $\chi^{2}$ Selection
Reconstructed $\pi^{0}\to\gamma\gamma$ in the final state using chi-square selection: 

$$\chi^2_{M_{\gamma\gamma}}=\frac{(M_{\gamma\gamma}-m_{\pi^{\text{0}}})^2}{\sigma^2_{M_{\gamma\gamma}}}, ~~~~~ \frac{\sigma_{M_{\gamma\gamma}}}{M_{\gamma\gamma}}=\frac{1}{2}\sqrt{\left(\frac{\sigma_{1}}{E_1}\right)^{2}+\left(\frac{\sigma_{2}}{E_2}\right)^{2}}$$

**Selection criteria:**
The pair of photons $\gamma_{1}$ and $\gamma_{2}$ is chosen such that the reconstructed mass $M_{\gamma\gamma}$ is closest to the mass constraint $$m_{\pi^{0}}=\sqrt{2E_{1}E_{2}\left(1-\text{cos}\theta_{12}\right)}.$$ 

The $\chi^{2}$-test is performed on event-by-event basis, and the energy-dependent relative error $\sigma_{M_{\gamma\gamma}}/M_{\gamma\gamma}$ is directly associated with the uncertainties of the photon energies $\sigma_{1}$, $\sigma_{2}$.\ Uncertainty contributions from both the measured angles and the correlations between $E_{1}$ and $E_{2}$ are considered negligible. The test is conducted for each photon pair, and the combination with the smallest chi-square value  $\chi^{2}_{M_{\gamma\gamma}}$ provides the best candidates for the $\pi^{0}$ decay photons. The efficiency of $\pi^{0}$ decay photon identification can be estimated by comparing the selected photons to the true MC information. 

**Related repository:** https://github.com/boaca926-beep/KLOE_REPO.git

# 📄 KLOE Analysis (BDT)

## 💡 Description
- Re-analyze KLOE analysis using BDT dataset.

- Reconstructed $\pi^{0}\to\gamma\gamma$ in the final state using BDT selection

**Related repository:** https://github.com/boaca926-beep/KLOE_BDT.git    
        
## 🚀  Data Preparation for BDT Selection

### 1. Input Raw Data Files
**Scripts**
```
script/listpath.sh # listing path of raw data root files stored as a text input file  
```

**Outputs**
```
path_chain/*
```

### 2. Create ROOT files (for both analysis and BDT training)
**Scripts**
```
run_bdt/Process.C (prompt, small samples)
run_bdt/script/input_bdt.sh (analysis, large samples)
```

**Outputs**
```
/home/bo/Desktop/sig.root (prompt)
/home/bo/Desktop/analysis_root_v6/input_bdt_TDATA_chain/input/sig.root (analysis)
```

### 3. Convert ROOT files to input for BDT training
**Scripts:**
```
run_bdt/get_bdt_sample.C        # prompt, small samples 
script/get_bdt_sample.sh        # analysis, large samples samples
```

**Outputs:**
```
dataset/sig_bdt.root
dataset/kloe_bdt.root
```

## 🚀 Obtain BDT Model (CUDA-boosted)
### Key Words
XGBoost, training, validation, test, model, ROC, AUC, confusion matrix

### Workflow Steps 
0. Setup Enviroment
    - Using UV pyproject.toml ```uv sync```
    
    - Using requirements.txt ```uv add -r requirements.txt```

    - Working Space ```/home/kloe/Desktop/KLOE_BDT/bdt```

1. Data preparation
    **Scripts:**
    ```
    # Spliting dataset to training, validate and test
    uv run main_initialize_kloe_opti.py \
        --input /home/kloe/Desktop/KLOE_BDT/dataset/kloe_bdt.root \
        --chunk-size 50000 \
        --output-dir /home/kloe/Desktop/KLOE_BDT/dataset_bdt
    ```

    **Outputs:**
    ```
    dataset_bdt/*
    ```

2. Inspect features - photons and photon pairs
    **Scripts:**
    ```
    # Inspecting photon features and features of all paried-photon combinations
    uv run main_inspect.py
    ```

    **Outputs:**
    ```
    /home/kloe/Desktop/KLOE_BDT/plots_inspect
    ```
    <!-- ![Kinematic Comparison: Photon vs Photon Pair Variables](plots_inspect/Kine_compr_TCOMB.png)
    *Figure 1: Comparison of kinematic variables between single photons and photon pairs from π⁰ decay*
    -->

    **Inspection Plots:**
    <div align="center">
    <img src="plots_inspect/Kine_compr_TCOMB.png" width="500" alt="Kinematic Variables"/>
    <br/>
    <em>Figure 1: Kinematic variables distribution after selection cuts</em>

    <br/><br/>

    <img src="plots_inspect/Photon_4-momentum_compr_TCOMB.png" width="500" alt="Photon Features"/>
    <br/>
    <em>Figure 2: 4-momentum features of photons in the final state</em>

    <br/><br/>

    <img src="plots_inspect/Pi0_compr_TCOMB.png" width="500" alt="Paired Photon Features"/>
    <br/>
    <em>Figure 3: Features of paired-photon combinations for π⁰ reconstruction</em>

    <br/><br/>

    <img src="plots_inspect/pos_pi0_mass_TCOMB.png" width="500" alt="Pi0 Mass"/>
    <br/>
    <em>Figure 4: Reconstructed π⁰ mass peak at nominal value (135 MeV/c²)</em>
    
    </div>
    

3. Training - hyperparameter tuning

4. Validate and test
    
    - Confusion matrix 
                   
    - ROC and AUC plots
    
    - Best threshold determination

5. Final testing

## 📊 Performance Metrics

- ROC curve and AUC score

- Confusion matrix

- Optimized threshold selection