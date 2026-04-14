# 📄 KLOE Analysis (BDT)

## 🚀  Data preparation for BDT selection (laptop)
- ### Workspace
    ```bash
    kloe@bo-ThinkPad-P16v-Gen-1:~/Desktop/KLOE_BDT$
    ```

- ### Create root files (for both analysis and BDT training)
    - Script
    ```bash
    run_bdt/Process.C (prompt, small samples)
    run_bdt/script/input_bdt.sh (analysis, large samples)
    ```

    - Output
    ```bash
    /home/bo/Desktop/analysis_root_v6/input_bdt_TDATA_chain/input/*.root (analysis)
    /home/kloe/Desktop/sig.root
    ```

- ### Convert root files to input for BDT training
    - Scripte
    ```bash
    run_bdt/get_bdt_sample.C (prompt, small samples)
    script/get_bdt_sample.sh (analysis, large smaples)
    ```

    - Output
    ```bash
    dataset/sig_bdt.root
    dataset/kloe_samples.root
    ```

## 🚀 Obtain BDT model (laptop, cuda gpu)

### Inspect features
- photons and photon pairs

### train the model
- hypertuning parameters

### validate and test
- confusion matrix, 
- ROC and AUC plots. 
- Best threshold