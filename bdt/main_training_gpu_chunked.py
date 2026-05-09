#!/usr/bin/env python3
# main_training_gpu_chunked.py
"""
GPU‑accelerated chunked training for arbitrarily large datasets.
Automatically converts existing .pkl files to .npy chunks if needed.
Requires XGBoost >= 3.0 and CuPy (for GPU chunk loading).
"""

import sys
import os
import glob
import gc
import time
import json
import shutil
import psutil
import multiprocessing
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import joblib
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from bayes_opti import baye_opti

from config import DATA_DIR, MODEL_DIR, patched_get_basescore

# ======================================================================
#  GPU and memory helpers
# ======================================================================
def has_gpu():
    try:
        import cupy as cp
        cp.cuda.runtime.getDeviceCount()
        return True
    except:
        return False

def check_memory_usage(threshold_gb=50):
    mem = psutil.virtual_memory()
    available_gb = mem.available / 1e9
    print(f"Available memory: {available_gb:.1f} GB")
    if available_gb < threshold_gb:
        warnings.warn(f"Only {available_gb:.1f}GB available. Large dataset may cause swapping.")
    return available_gb

# ======================================================================
#  Chunk conversion from pickle files
# ======================================================================
def convert_pickles_to_chunks(br_type, chunk_size=100_000):
    """Convert X_train, y_train, X_val, y_val .pkl files into chunked .npy files."""
    print("\n🔄 Converting pickled data to chunked .npy files ...")
    # Load pickles
    try:
        X_train = joblib.load(os.path.join(DATA_DIR, f'X_train_{br_type}.pkl'))
        y_train = joblib.load(os.path.join(DATA_DIR, f'y_train_{br_type}.pkl'))
        X_val   = joblib.load(os.path.join(DATA_DIR, f'X_val_{br_type}.pkl'))
        y_val   = joblib.load(os.path.join(DATA_DIR, f'y_val_{br_type}.pkl'))
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing pickle file: {e}. Run main_initialize_kloe_opti.py first.")

    # Convert to numpy arrays
    if hasattr(X_train, 'values'):
        X_train = X_train.values
    if hasattr(y_train, 'values'):
        y_train = y_train.values
    if hasattr(X_val, 'values'):
        X_val = X_val.values
    if hasattr(y_val, 'values'):
        y_val = y_val.values
    y_train = y_train.ravel()
    y_val = y_val.ravel()

    # Save feature names
    try:
        X_train_df = joblib.load(os.path.join(DATA_DIR, f'X_train_{br_type}.pkl'))
        if hasattr(X_train_df, 'columns'):
            joblib.dump(X_train_df.columns.tolist(),
                        os.path.join(DATA_DIR, f'training_cols_{br_type}.pkl'))
    except:
        pass

    # Helper to write chunks
    def save_chunks(data, prefix, chunk_size):
        n = len(data)
        n_chunks = (n + chunk_size - 1) // chunk_size
        for i in range(n_chunks):
            start = i * chunk_size
            end = min(start + chunk_size, n)
            chunk = data[start:end].astype(np.float32 if 'X' in prefix else np.uint8)
            out_path = os.path.join(DATA_DIR, f'{prefix}_{br_type}_chunk_{i}.npy')
            np.save(out_path, chunk)
            print(f"  Saved {out_path} ({len(chunk)} samples)")

    # Remove any existing old .npy chunks to avoid mixing
    for pattern in ['X_train', 'y_train', 'X_val', 'y_val']:
        for f in glob.glob(os.path.join(DATA_DIR, f'{pattern}_{br_type}_chunk_*.npy')):
            os.remove(f)

    print("Writing training chunks...")
    save_chunks(X_train, 'X_train', chunk_size)
    save_chunks(y_train, 'y_train', chunk_size)

    print("Writing validation chunks...")
    save_chunks(X_val, 'X_val', chunk_size)
    save_chunks(y_val, 'y_val', chunk_size)

    print("✅ Conversion finished.\n")

def ensure_chunks(br_type, chunk_size=100_000):
    """Check if chunk files exist; if not, run conversion."""
    mandatory = [f'X_train_{br_type}_chunk_0.npy', f'y_train_{br_type}_chunk_0.npy',
                 f'X_val_{br_type}_chunk_0.npy', f'y_val_{br_type}_chunk_0.npy']
    for f in mandatory:
        if not os.path.exists(os.path.join(DATA_DIR, f)):
            convert_pickles_to_chunks(br_type, chunk_size)
            return

# ======================================================================
#  Subset loader for Bayesian optimisation
# ======================================================================
def load_dataset_subset(br_type, sample_fraction=0.1):
    """Load a small subset of chunked data for Bayesian optimisation."""
    input_dir = DATA_DIR
    X_chunks = sorted(glob.glob(os.path.join(input_dir, f'X_train_{br_type}_chunk_*.npy')))
    y_chunks = sorted(glob.glob(os.path.join(input_dir, f'y_train_{br_type}_chunk_*.npy')))
    n_chunks = max(1, int(len(X_chunks) * sample_fraction))
    X_list, y_list = [], []
    total_samples = 0
    for i in range(n_chunks):
        X = np.load(X_chunks[i])
        y = np.load(y_chunks[i]).ravel()
        X_list.append(X)
        y_list.append(y)
        total_samples += X.shape[0]
        if total_samples > 500000:
            break
    X_sample = np.vstack(X_list) if len(X_list) > 1 else X_list[0]
    y_sample = np.hstack(y_list) if len(y_list) > 1 else y_list[0]
    # feature names
    training_cols_path = os.path.join(input_dir, f'training_cols_{br_type}.pkl')
    if os.path.exists(training_cols_path):
        training_cols = joblib.load(training_cols_path)
        X_sample = pd.DataFrame(X_sample, columns=training_cols)
    else:
        training_cols = [f'feature_{i}' for i in range(X_sample.shape[1])]
    print(f"Loaded {len(X_sample)} samples for Bayesian optimization")
    return X_sample, y_sample, training_cols

# ======================================================================
#  Custom iterator for GPU chunked loading
# ======================================================================
class GPUChunkedDataIter(xgb.DataIter):
    def __init__(self, X_files, y_files, cache_prefix="./gpu_cache"):
        super().__init__(cache_prefix=cache_prefix)
        self.X_files = X_files
        self.y_files = y_files
        self.it = 0

    def _load_chunk(self, idx):
        if idx >= len(self.X_files):
            return None, None
        import cupy as cp
        X = cp.load(self.X_files[idx])   # load directly as GPU array
        y = cp.load(self.y_files[idx]).ravel()
        return X, y

    def next(self, input_data):
        if self.it >= len(self.X_files):
            return False
        X, y = self._load_chunk(self.it)
        input_data(data=X, label=y)
        self.it += 1
        return True

    def reset(self):
        self.it = 0

# ======================================================================
#  Main training
# ======================================================================
if __name__ == '__main__':
    print("GPU + chunked training for large datasets")
    print("=" * 60)
    check_memory_usage()
    gpu_available = has_gpu()
    if gpu_available:
        print("✓ GPU detected – will use external memory GPU training.")
        import cupy as cp
        cp.cuda.set_allocator(cp.cuda.MemoryPool().malloc)
    else:
        print("⚠ No GPU found – falling back to CPU chunked training.")
        # (ExtMemQuantileDMatrix also works on CPU)

    # --- Prepare paths and ensure chunks exist ---
    model_dir = MODEL_DIR
    os.makedirs(model_dir, exist_ok=True)
    input_dir = DATA_DIR

    br_type = 'TCOMB'
    phys_map = joblib.load(os.path.join(input_dir, 'phys_map.pkl'))
    info = phys_map[br_type]
    print(f"Training dataset: {br_type} - {info['br_title']}")

    # Automatically create chunks if missing
    ensure_chunks(br_type, chunk_size=100_000)

    # --- Collect chunk file lists ---
    X_chunks = sorted(glob.glob(os.path.join(input_dir, f'X_train_{br_type}_chunk_*.npy')))
    y_chunks = sorted(glob.glob(os.path.join(input_dir, f'y_train_{br_type}_chunk_*.npy')))
    X_val_chunks = sorted(glob.glob(os.path.join(input_dir, f'X_val_{br_type}_chunk_*.npy')))
    y_val_chunks = sorted(glob.glob(os.path.join(input_dir, f'y_val_{br_type}_chunk_*.npy')))

    print(f"Found {len(X_chunks)} training chunks, {len(X_val_chunks)} validation chunks.")

    # --- Bayesian optimisation on a subset ---
    print("\nSTEP 1: Bayesian optimisation on 10% subset")
    X_subset, y_subset, feature_names = load_dataset_subset(br_type, sample_fraction=0.1)
    if gpu_available:
        X_subset = cp.array(X_subset)
        y_subset = cp.array(y_subset)
    X_subset_np = cp.asnumpy(X_subset) if gpu_available else X_subset
    y_subset_np = cp.asnumpy(y_subset) if gpu_available else y_subset
    optimized_params = baye_opti(X_subset_np, y_subset_np)
    print("Optimised hyperparameters:", optimized_params)

    # --- Training parameters ---
    params = {
        'tree_method': 'hist',          # required for GPU external memory
        'device': 'cuda' if gpu_available else 'cpu',
        'nthread': multiprocessing.cpu_count(),
        'max_bin': 512,
        'eval_metric': ['auc', 'logloss'],
        'early_stopping_rounds': 50,
        'verbosity': 1,
        'max_depth': optimized_params.get('max_depth', 10),
        'learning_rate': optimized_params.get('learning_rate', 0.1),
        'subsample': optimized_params.get('subsample', 0.8),
        'colsample_bytree': optimized_params.get('colsample_bytree', 0.8),
        'min_child_weight': optimized_params.get('min_child_weight', 1),
        'gamma': optimized_params.get('gamma', 0),
        'reg_alpha': optimized_params.get('reg_alpha', 0),
        'reg_lambda': optimized_params.get('reg_lambda', 1),
    }

    print("\nSTEP 2: Creating external memory DMatrix with GPU streaming")
    train_iter = GPUChunkedDataIter(X_chunks, y_chunks, cache_prefix="./cache_train")
    dtrain = xgb.ExtMemQuantileDMatrix(train_iter)

    if X_val_chunks:
        val_iter = GPUChunkedDataIter(X_val_chunks, y_val_chunks, cache_prefix="./cache_val")
        dval = xgb.ExtMemQuantileDMatrix(val_iter)
        evals = [(dtrain, 'train'), (dval, 'val')]
    else:
        evals = [(dtrain, 'train')]

    print("\nSTEP 3: Starting GPU‑accelerated external memory training")
    start_time = time.time()
    evals_result = {}
    booster = xgb.train(
        params,
        dtrain,
        num_boost_round=1000,
        evals=evals,
        evals_result=evals_result,
        early_stopping_rounds=params['early_stopping_rounds'],
        verbose_eval=50
    )
    training_time = time.time() - start_time

    # Extract best AUC
    if 'val' in evals_result and 'auc' in evals_result['val']:
        best_auc = evals_result['val']['auc'][booster.best_iteration]
    else:
        best_auc = evals_result['train']['auc'][booster.best_iteration]

    # --- Save model and metrics ---
    n_features = len(feature_names)
    model_wrapper = xgb.XGBClassifier(**params)
    model_wrapper._Booster = booster
    joblib.dump(model_wrapper, f'{model_dir}/pi0_classifier_model_{br_type}.pkl', compress=3)
    booster.save_model(f'{model_dir}/pi0_classifier_model_{br_type}.json')

    metrics = {
        'auc': float(best_auc),
        'best_iteration': booster.best_iteration,
        'best_score': booster.best_score,
        'params': params,
        'n_features': n_features,
        'feature_names': feature_names,
        'training_time_minutes': training_time / 60,
        'gpu_enabled': gpu_available,
        'external_memory': True,
    }
    with open(f'{model_dir}/metrics_{br_type}.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"✅ AUC on validation: {best_auc:.4f}")
    print(f"Training finished in {training_time/60:.1f} minutes")

    # --- Save to ROOT format for TMVA ---
    try:
        import ROOT
        try:
            import ROOT._pythonization._tmva._tree_inference as tree_inference
            tree_inference.get_basescore = patched_get_basescore
        except:
            pass
        ROOT.TMVA.Experimental.SaveXGBoost(
            booster,
            "BDT_pi0",
            f"{model_dir}/bdt_pi0_{br_type}.root",
            num_inputs=n_features
        )
        print(f"✓ ROOT model saved to {model_dir}/bdt_pi0_{br_type}.root")
    except Exception as e:
        print(f"⚠ Could not save ROOT model: {e}")

    # Cleanup
    del dtrain
    if 'val_iter' in locals(): del val_iter
    gc.collect()
    print("\n✓ GPU chunked training completed.")