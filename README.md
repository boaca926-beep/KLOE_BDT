# 📄 KLOE Analysis (Standard)
## 💡 Description
- ROOT framework based analysis.

- Study of $e^{+}e^{-}\to\pi^{+}\pi^{-}\pi^{0}\gamma\to\pi^{+}\pi^{-}\gamma\gamma\gamma$ ISR process.

- Reconstructed $\pi^{0}\to\gamma\gamma$ in the final state using chi-square selection. $$\chi^2_{M_{\gamma\gamma}}=\frac{(M_{\gamma\gamma}-m_{\pi^{\text{0}}})^2}{\sigma^2_{M_{\gamma\gamma}}}, ~~~~~ \frac{\sigma_{M_{\gamma\gamma}}}{M_{\gamma\gamma}}=\frac{1}{2}\sqrt{\left(\frac{\sigma_{1}}{E_1}\right)^{2}+\left(\frac{\sigma_{2}}{E_2}\right)^{2}}$$

    The pair of photons $\gamma_{1}$ and $\gamma_{2}$ is chosen such that the reconstructed mass $M_{\gamma\gamma}$ is closest to the mass constraint 
    $$m_{\pi^{0}}=\sqrt{2E_{1}E_{2}\left(1-\text{cos}\theta_{12}\right)}.$$ 

    The $\chi^{2}$-test is performed on event-by-event basis, and the energy-dependent relative error $\sigma_{M_{\gamma\gamma}}/M_{\gamma\gamma}$ is directly associated with the uncertainties of the photon energies $\sigma_{1}$, $\sigma_{2}$.\ Uncertainty contributions from both the measured angles and the correlations between $E_{1}$ and $E_{2}$ are considered negligible. The test is conducted for each photon pair, and the combination with the smallest chi-square value  $\chi^{2}_{M_{\gamma\gamma}}$ provides the best candidates for the $\pi^{0}$ decay photons. The efficiency of $\pi^{0}$ decay photon identification can be estimated by comparing the selected photons to the true MC information. 

- https://github.com/boaca926-beep/KLOE_REPO.git

# 📄 KLOE Analysis (BDT)
## 💡 Description
- Re-analyze KLOE analysis using BDT dataset.

- Reconstructed $\pi^{0}\to\gamma\gamma$ in the final state using BDT selection. 

- https://github.com/boaca926-beep/KLOE_BDT.git    
        
## 🚀  Data preparation for BDT selection

- ### Create root files (for both analysis and BDT training)
    - Script
    ```bash
    run_bdt/Process.C (prompt, small samples)
    run_bdt/script/input_bdt.sh (analysis, large samples)
    ```

    - Output
    ```bash
    /home/bo/Desktop/sig.root (promt)
    /home/bo/Desktop/analysis_root_v6/input_bdt_TDATA_chain/input/sig.root (analysis)
    ```

- ### Convert root files to input for BDT training
    - Script
    ```bash
    run_bdt/get_bdt_sample.C (prompt, small samples)
    script/get_bdt_sample.sh (analysis, large smaples)
    ```

    - Output
    ```bash
    dataset/sig_bdt.root
    dataset/kloe_bdt.root
    ```

## 🚀 Obtain BDT model (cuda-boosted)
- ### Key words
    - XGBoost, training, valication, test, model, ROC, AUC, confusion matrix

- ### Inspect features
    - photons and photon pairs

- ### Training
    - hypertuning parameters

- ### Validate and test
    - confusion matrix, 
    
    - ROC and AUC plots. 
    
    - Best threshold

- ### Test
