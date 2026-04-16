# Training script
import uproot
import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

import numpy as np
import pandas as pd
from bayes_opti import baye_opti
import joblib
import xgboost as xgb
import ROOT
import json
import multiprocessing
import psutil
import time
import gc
import glob

from config import (
    DATA_DIR,
    MODEL_DIR,
    patched_get_basescore
)

'''
Chunked approach is much slower than the in-memory approach, but barely works with arbitrarily large dataset, and insufficient RAM size 

Version	Time	CPU	Best For
In-memory (float32)	6 seconds	7.8%	Small datasets (< 2GB)
Chunked	72 seconds	6.5%	Large datasets (> 16GB)

Dataset Size    RAM Required    Works?
─────────────────────────────────────────
10 GB           20 GB RAM       ✅ Yes (in-memory)
50 GB           100 GB RAM      ❌ No (would crash)
50 GB           16 GB RAM       ✅ Yes (chunked) ← This is the point!
500 GB          16 GB RAM       ✅ Yes (chunked)
1 TB            16 GB RAM       ✅ Yes (chunked)
'''

def load_dataset(br_type):
    """
    Load train, val dataset - MODIFIED to handle chunked data
    """
    global input_data_dir
    input_data_dir_local = input_data_dir
    
    # Check if chunked data exists
    chunk_pattern = os.path.join(input_data_dir_local, f'X_train_{br_type}_chunk_*.npy')
    has_chunks = len(glob.glob(chunk_pattern)) > 0
    
    if has_chunks:
        print("Loading chunked data (external memory mode)...")
        return None, None, None, None
    else:
        # Original loading for smaller datasets
        X_train = joblib.load(os.path.join(input_data_dir_local, f'X_train_{br_type}.pkl'))
        y_train = joblib.load(os.path.join(input_data_dir_local, f'y_train_{br_type}.pkl'))
        X_val = joblib.load(os.path.join(input_data_dir_local, f'X_val_{br_type}.pkl'))
        y_val = joblib.load(os.path.join(input_data_dir_local, f'y_val_{br_type}.pkl'))
        return X_train, y_train, X_val, y_val


def load_dataset_subset(br_type, sample_fraction=0.1):
    """
    Load subset of chunked data for Bayesian optimization
    """
    global input_data_dir
    input_data_dir_local = input_data_dir
    
    X_chunks = sorted(glob.glob(os.path.join(input_data_dir_local, f'X_train_{br_type}_chunk_*.npy')))
    y_chunks = sorted(glob.glob(os.path.join(input_data_dir_local, f'y_train_{br_type}_chunk_*.npy')))
    
    n_chunks_to_load = max(1, int(len(X_chunks) * sample_fraction))
    X_list = []
    y_list = []
    
    for i in range(n_chunks_to_load):
        X_chunk = np.load(X_chunks[i])
        y_chunk = np.load(y_chunks[i])
        X_list.append(X_chunk)
        y_list.append(y_chunk.ravel())
        
        # Limit to 500k samples to prevent memory issues
        if sum(x.shape[0] for x in X_list) > 500000:
            break
    
    X_sample = np.vstack(X_list) if len(X_list) > 1 else X_list[0]
    y_sample = np.hstack(y_list) if len(y_list) > 1 else y_list[0]
    
    # Load feature names if available
    training_cols_path = os.path.join(input_data_dir_local, f'training_cols_{br_type}.pkl')
    if os.path.exists(training_cols_path):
        training_cols = joblib.load(training_cols_path)
        X_sample = pd.DataFrame(X_sample, columns=training_cols)
    
    print(f"Loaded {len(X_sample)} samples for Bayesian optimization")
    return X_sample, y_sample, None, None


if __name__ == '__main__':

    print(f"Train individual signal physical channels ...")
    model_dir = MODEL_DIR #"./models" 

    import shutil
    if os.path.exists(model_dir):
        shutil.rmtree(model_dir)

    os.makedirs(model_dir, exist_ok=True)

    input_data_dir = DATA_DIR
    phys_map = joblib.load(os.path.join(input_data_dir, f'phys_map.pkl'))
    print(phys_map)

    br_type = 'TCOMB'
    info = phys_map[br_type]
    info_title = info['br_title']
    info_category = info['category']
    print(f"Inspecting dataset of data_type: {br_type}; title: {info_title}; category: {info_category}")  
    
    # Check if we have chunked data
    chunk_pattern = os.path.join(input_data_dir, f'X_train_{br_type}_chunk_*.npy')
    has_chunks = len(glob.glob(chunk_pattern)) > 0
    
    print("\n" + "="*60)
    print("STEP 1: BAYESIAN OPTIMIZATION")
    print("="*60)

    # MODIFIED: Use subset for Bayesian optimization if data is large
    if has_chunks:
        print("Large dataset detected - using 10% subset for Bayesian optimization...")
        X_train_subset, y_train_subset, _, _ = load_dataset_subset(br_type, sample_fraction=0.1)
        optimized_params = baye_opti(X_train_subset, y_train_subset)
        del X_train_subset, y_train_subset
        gc.collect()
        
        # Load feature columns from saved metadata
        training_cols_path = os.path.join(input_data_dir, f'training_cols_{br_type}.pkl')
        if os.path.exists(training_cols_path):
            training_cols = joblib.load(training_cols_path)
        else:
            # If no metadata, load from first chunk
            X_first = np.load(glob.glob(os.path.join(input_data_dir, f'X_train_{br_type}_chunk_*.npy'))[0])
            training_cols = [f'feature_{i}' for i in range(X_first.shape[1])]
    else:
        # Original path for smaller datasets
        X_train, y_train, X_val, y_val = load_dataset(br_type)
        print(f"X_train columns: {X_train.columns}")
        print(f"Number of features in the training: {len(X_train.columns.tolist())}")
        optimized_params = baye_opti(X_train, y_train)
        training_cols = [col for col in X_train.columns if col != 'is_pi0']

    # Set fixed parameters
    params = {
         'nthread': -1,
         'tree_method': 'hist',
         #'early_stopping_rounds': 50,
         'eval_metric': ['auc', 'error'],
         'verbosity': 1,
         'max_depth': optimized_params.get('max_depth', 10),
         'learning_rate': optimized_params.get('learning_rate', 0.1),
         'subsample': optimized_params.get('subsample', 0.8),
         'colsample_bytree': optimized_params.get('colsample_bytree', 0.8),
         'min_child_weight': optimized_params.get('min_child_weight', 1),
         'gamma': optimized_params.get('gamma', 0),
         'reg_alpha': optimized_params.get('reg_alpha', 0),
         'reg_lambda': optimized_params.get('reg_lambda', 1)
         }
    
    print("\n" + "="*50)
    print("THREAD CONFIGURATION CHECK")
    print("="*50)
    print(f"CPU cores available: {multiprocessing.cpu_count()}")
    print(f"XGBoost nthread setting: {params.get('nthread')}")
    print(f"XGBoost version: {xgb.__version__}")

    print(f"\n{'='*60}")
    print(f"STEP 2. TRAIN BEST MODEL")
    print(f"{'='*60}")    

    # MODIFIED: Handle training differently based on data size
    if has_chunks:
        # EXTERNAL MEMORY TRAINING FOR LARGE DATASET
        print("\nUsing external memory training for large dataset...")
        
        # Get all chunk files
        X_chunks = sorted(glob.glob(os.path.join(input_data_dir, f'X_train_{br_type}_chunk_*.npy')))
        y_chunks = sorted(glob.glob(os.path.join(input_data_dir, f'y_train_{br_type}_chunk_*.npy')))
        X_val_chunks = sorted(glob.glob(os.path.join(input_data_dir, f'X_val_{br_type}_chunk_*.npy')))
        y_val_chunks = sorted(glob.glob(os.path.join(input_data_dir, f'y_val_{br_type}_chunk_*.npy')))
        
        # FIX: Create a single pattern string for external memory
        # XGBoost expects a string with wildcards, not a list of strings
        # The pattern should match all chunk files
        train_pattern = os.path.join(input_data_dir, f'X_train_{br_type}_chunk_*.npy')
        
        # For external memory with separate label files, we need to use the libsvm format
        # or create a temporary directory with the correct structure
        # Simpler approach: Create a single file list using the correct URI format for each chunk
        # but XGBoost expects a single string with # symbol for multiple files
        
        # ALTERNATIVE: Use DMatrix with custom iterator (works with XGBoost 1.7+)
        # Let's use the DataIter approach which is more reliable
        
        from xgboost import DataIter
        
        class ChunkedDataIter(DataIter):
            def __init__(self, X_files, y_files):
                self.X_files = X_files
                self.y_files = y_files
                self.it = 0
                super().__init__()
            
            def next(self, input_data):
                if self.it >= len(self.X_files):
                    return False
                
                # Load chunk
                X = np.load(self.X_files[self.it])
                y = np.load(self.y_files[self.it]).ravel()
                
                # Pass to XGBoost
                input_data(data=X, label=y)
                
                self.it += 1
                return True
            
            def reset(self):
                self.it = 0
        
        # Create iterator
        train_iter = ChunkedDataIter(X_chunks, y_chunks)
        dtrain = xgb.DMatrix(train_iter)
        
        # For validation
        if X_val_chunks:
            val_iter = ChunkedDataIter(X_val_chunks, y_val_chunks)
            dval = xgb.DMatrix(val_iter)
            evals = [(dtrain, 'train'), (dval, 'val')]
        else:
            evals = [(dtrain, 'train')]
        
        print(f"Training chunks: {len(X_chunks)}")
        print(f"Validation chunks: {len(X_val_chunks) if X_val_chunks else 0}")
        
        # Train with external memory
        start_time = time.time()
        evals_result = {}
        booster = xgb.train(
            params,
            dtrain,
            num_boost_round=1000,
            evals=evals,
            evals_result=evals_result,
            early_stopping_rounds=50,
            verbose_eval=50
        )
        training_time = time.time() - start_time
        
        # Create wrapper for compatibility
        model = xgb.XGBClassifier(**params)
        model._Booster = booster
        
        # Get best AUC
        if 'val' in evals_result and 'auc' in evals_result['val']:
            best_auc = evals_result['val']['auc'][booster.best_iteration]
        else:
            best_auc = evals_result['train']['auc'][booster.best_iteration]
        
        # Get number of features
        n_features = len(training_cols)
        
        # Save model
        booster.save_model(f'{model_dir}/pi0_classifier_model_{br_type}.json')
        joblib.dump(model, f'{model_dir}/pi0_classifier_model_{br_type}.pkl')
        
        print(f"\nTraining completed in {training_time/60:.1f} minutes")
        
    else:
        # ORIGINAL IN-MEMORY TRAINING
        X_train, y_train, X_val, y_val = load_dataset(br_type)
        
        X_train_np = X_train[training_cols].to_numpy()
        X_val_np = X_val[training_cols].to_numpy()
        y_train_np = np.asarray(y_train).ravel()
        y_val_np = np.asarray(y_val).ravel()

        model = xgb.XGBClassifier(
            **params,
            early_stopping_rounds=50
            ) 
        print(f"Model will use: {model.get_params().get('nthread', 'default')} threads")
        
        print("\nTraining started - monitoring CPU...")
        start_time = time.time()

        model.fit(X_train_np, y_train_np,
                  eval_set=[(X_train_np, y_train_np), (X_val_np, y_val_np)],
                  verbose=50)
        
        training_time = time.time() - start_time
        
        from sklearn.metrics import roc_auc_score
        y_pred = model.predict_proba(X_val_np)[:, 1]
        best_auc = roc_auc_score(y_val_np, y_pred)
        
        joblib.dump(model, f'{model_dir}/pi0_classifier_model_{br_type}.pkl')
        n_features = len(training_cols)
        
        del X_train_np, X_val_np, y_train_np, y_val_np
        
        booster = model.get_booster()

    # CPU usage monitoring
    cpu_percent = psutil.cpu_percent(interval=1, percpu=True)
    print(f"\nCPU usage per core during training: {cpu_percent}")
    print(f"Average CPU usage: {sum(cpu_percent)/len(cpu_percent):.1f}%")
    print(f"Peak threads used: {psutil.Process().num_threads()}")
    print("="*50)

    print(f"\n{'='*60}")
    print(f"STEP 3: SAVE METRICS")
    print(f"{'='*60}") 
    
    # Save metrics
    metrics = {
         'auc': float(best_auc),
         'best_iteration': booster.best_iteration if has_chunks else model.best_iteration,
         'best_score': booster.best_score if has_chunks else model.best_score,
         'params': params,
         'n_features': n_features,
         'training_time_minutes': training_time / 60,
         'external_memory': has_chunks
         }
            
    with open(f'{model_dir}/metrics_{br_type}.json', 'w') as f:
        json.dump(metrics, f, indent=2)
        print(f"Metrics saved to {model_dir}/metrics_{br_type}.json")
        print(f"AUC: {metrics['auc']:.4f}")

    # Save to ROOT format
    import ROOT._pythonization._tmva._tree_inference as tree_inference
    tree_inference.get_basescore = patched_get_basescore

    # Get the booster for ROOT saving
    if not has_chunks:
        booster = model.get_booster()
    
    try:
        ROOT.TMVA.Experimental.SaveXGBoost(
            model,  # ← This is XGBClassifier, which has .objective attribute
            "BDT_pi0", 
            f"{model_dir}/bdt_pi0_{br_type}.root", 
            num_inputs=n_features
        )
        print(f"✓ Model saved to {model_dir}/bdt_pi0_{br_type}.root")
    except Exception as e:
        print(f"✗ Failed to save ROOT model: {e}")

    # Memory cleanup
    gc.collect()
            
    print(f"\n✓ Training completed for {br_type}")