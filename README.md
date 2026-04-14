# KLOE Analysis (BDT)

## Data preparation for BDT selection (laptop)
- Location: kloe@bo-ThinkPad-P16v-Gen-1:~/Desktop/KLOE_BDT$

### Create root files (for both analysis and BDT training)
- run_bdt/Process.C (prompt, small samples)
- run_bdt/script/input_bdt.sh (analysis, large samples)

OUTPUT
```bash
- /home/bo/Desktop/analysis_root_v6/input_bdt_TDATA_chain/input/*.root (analysis)
- /home/kloe/Desktop/sig.root
```

### Convert OUTPUT/*root to sample root files (for BDT training)
```bash
- run_bdt/get_bdt_sample.C (prompt, small samples)
output: /home/kloe/Desktop/KLOE_BDT/dataset/sig_bdt.root
```

```bash
- script/get_bdt_sample.sh (analysis, large smaples)
output: /home/kloe/Desktop/KLOE_BDT/dataset/kloe_samples.root
```

## Obtain BDT model (laptop, cuda gpu)

### Inspect features: photons and photon pairs
- train the model: hypertuning parameters
- validate and test: confusion matrix, ROC and AUC plots. Best threshold