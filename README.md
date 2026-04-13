# KLOE Analysis (BDT)

## Data preparation for BDT selection (laptop)
- Location: kloe@bo-ThinkPad-P16v-Gen-1:~/Desktop/KLOE_BDT$

### Create root files (for both analysis and BDT training)
- run_bdt/Process.C (for small samples)
- run_bdt/script/input_bdt.sh (for intermidate or large samples)
```bash
OUTPUT: /home/bo/Desktop/analysis_root_v6/input_bdt_TDATA_chain/input/*.root
```

### Convert OUTPUT/*root to sample root files (for BDT training)
- script/get_sample.sh

## Obtain BDT model (laptop, cuda gpu)

### Inspect features: photons and photon pairs
- train the model: hypertuning parameters
- validate and test: confusion matrix, ROC and AUC plots. Best threshold