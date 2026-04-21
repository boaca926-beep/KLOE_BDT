# Training script - Memory optimized version (GPU ADAPTED)
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
import warnings

from config import (
    DATA_DIR, MODEL_DIR,
    patched_get_basescore
)


def check_memory_usage(threshold_gb=50):
    """Check if we have enough memory before loading"""
    mem = psutil.virtual_memory()
    available_gb = mem.available / 1e9
    print(f"Available memory: {available_gb:.1f} GB")
    
    if available_gb < threshold_gb:
        warnings.warn(f"Only {available_gb:.1f}GB available. Loading {threshold_gb}GB dataset may cause swapping.")
    return available_gb

# ========== ADDED: GPU check function ==========
def check_gpu():
    """Check if GPU is available for XGBoost"""
    try:
        # Test if GPU is available
        import subprocess
        result = subprocess.run(['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            gpu_name = result.stdout.strip()
            print(f"✓ GPU detected: {gpu_name}")
            
            # Check XGBoost GPU support
            if xgb.__version__ >= '1.6.0':
                print(f"✓ XGBoost {xgb.__version__} supports GPU")
                return True
            else:
                print(f"⚠ XGBoost {xgb.__version__} - upgrade to 1.6+ for GPU")
                return False
    except:
        pass
    
    print("⚠ No GPU detected or nvidia-smi not found - using CPU")
    return False
# ==============================================

def load_dataset_optimized(br_type):
    """
    Load train, val dataset with memory optimization
    """
    input_data_dir_local = input_data_dir
    
    # Load with memory optimization
    print("Loading X_train...")
    X_train = joblib.load(os.path.join(input_data_dir_local, f'X_train_{br_type}.pkl'))
    
    print("Loading y_train...")
    y_train = joblib.load(os.path.join(input_data_dir_local, f'y_train_{br_type}.pkl'))
    
    print("Loading X_val...")
    X_val = joblib.load(os.path.join(input_data_dir_local, f'X_val_{br_type}.pkl'))
    
    print("Loading y_val...")
    y_val = joblib.load(os.path.join(input_data_dir_local, f'y_val_{br_type}.pkl'))
    
    # Convert to more memory-efficient types if needed
    for col in X_train.select_dtypes(include=['float64']).columns:
        X_train[col] = X_train[col].astype('float32')
    
    for col in X_val.select_dtypes(include=['float64']).columns:
        X_val[col] = X_val[col].astype('float32')
    
    return X_train, y_train, X_val, y_val

def free_memory():
    """Force garbage collection and check memory"""
    gc.collect()
    mem = psutil.virtual_memory()
    print(f"Memory after cleanup: {mem.used/1e9:.1f}GB / {mem.total/1e9:.1f}GB ({mem.percent}%)")
    return mem

if __name__ == '__main__':

    print(f"Train individual signal physical channels ...")
    
    # Check memory before starting
    print("\n" + "="*50)
    print("MEMORY CHECK")
    print("="*50)
    initial_mem = check_memory_usage(threshold_gb=50)
    
    # ========== ADDED: Check GPU ==========
    print("\n" + "="*50)
    print("GPU CHECK")
    print("="*50)
    has_gpu = check_gpu()
    # ======================================
    
    # Create output directory
    model_dir = MODEL_DIR #"./models" 

    import shutil
    if os.path.exists(model_dir):
        shutil.rmtree(model_dir)

    os.makedirs(model_dir, exist_ok=True)

    input_data_dir = DATA_DIR
    phys_map = joblib.load(os.path.join(input_data_dir, f'phys_map.pkl'))
    print(phys_map)

    br_type = 'TCOMB' #'TETAGAM' 'TCOMB'
    info = phys_map[br_type]
    info_title = info['br_title']
    info_category = info['category']
    print(f"Inspecting dataset of data_type: {br_type}; title: {info_title}; category: {info_category}")  
    
    print("\n" + "="*60)
    print("STEP 1: BAYESIAN OPTIMIZATION")
    print("="*60)

    # Load dataset with memory optimization
    X_train, y_train, X_val, y_val = load_dataset_optimized(br_type)
    
    # Check memory after loading
    mem_after_load = free_memory()
    
    print(f"X_train shape: {X_train.shape}, memory: {X_train.memory_usage(deep=True).sum()/1e9:.2f}GB")
    print(f"X_train columns: {X_train.columns.tolist()[:10]}...")  # Show first 10
    print(f"Number of features: {len(X_train.columns.tolist())}")

    # Bayesian optimization (this will use some memory)
    print("\nRunning Bayesian optimization...")
    optimized_params = baye_opti(X_train, y_train)
    
    # Free memory from optimization temporary objects
    free_memory()

    # ========== MODIFIED: Parameters with GPU support ==========
    params = {
        # GPU settings optimized for 35 GPU
        #'tree_method': 'gpu_hist' if has_gpu else 'hist',  # GPU if available
        'tree_method': 'hist',  # 'hist' automatically uses GPU when device='cuda'
        'device': 'cuda' if has_gpu else 'cpu',            # Explicit device
        #'gpu_id': 0, # XGBoost Version < 3.0
        'nthread': 12,                      # -1: Use all available threads, don't oversubscribe on low-power GPU
        'max_bin': 512,                    # Reduce meomory footrpint
        #'gpu_hist_use_bfloat16': True,     # Use half-precision for efficiency, XGBoost auto-optimizes this

        # Training settings
        'early_stopping_rounds': 50,        # Early stopping
        'eval_metric': ['auc', 'error'],    # Evaluation metric
        'verbosity': 1,                     # Show progress
        'max_depth': optimized_params.get('max_depth', 10),
        'learning_rate': optimized_params.get('learning_rate', 0.1),
        'subsample': optimized_params.get('subsample', 0.8),
        'colsample_bytree': optimized_params.get('colsample_bytree', 0.8),
        'min_child_weight': optimized_params.get('min_child_weight', 1),
        'gamma': optimized_params.get('gamma', 0),
        'reg_alpha': optimized_params.get('reg_alpha', 0),
        'reg_lambda': optimized_params.get('reg_lambda', 1)
        }
    # ============================================================
    
    print("\n" + "="*50)
    print("CONFIGURATION CHECK")
    print("="*50)
    print(f"CPU cores available: {multiprocessing.cpu_count()}")
    print(f"XGBoost version: {xgb.__version__}")
    print(f"Tree method: {params.get('tree_method')}")      # ADDED
    print(f"Device: {params.get('device')}")                # ADDED
    print(f"nthread setting: {params.get('nthread')}")

    print(f"\n{'='*60}")
    print(f"STEP 2. TRAIN BEST MODEL")
    print(f"{'='*60}")    

    # Define training columns
    training_cols = [col for col in X_train.columns if col != 'is_pi0']

    # Convert to numpy arrays with float32 to save memory
    print("\nConverting to numpy arrays (float32)...")
    X_train_np = X_train[training_cols].to_numpy().astype(np.float32)
    X_val_np = X_val[training_cols].to_numpy().astype(np.float32)
    y_train_np = np.asarray(y_train).ravel().astype(np.float32)
    y_val_np = np.asarray(y_val).ravel().astype(np.float32)
    
    # Free pandas DataFrames to save memory
    del X_train, y_train, X_val, y_val
    free_memory()
    
    print(f"X_train_np shape: {X_train_np.shape}, dtype: {X_train_np.dtype}")
    print(f"Memory usage: {X_train_np.nbytes/1e9:.2f}GB (training), {X_val_np.nbytes/1e9:.2f}GB (validation)")

    # Create model with GPU settings
    model = xgb.XGBClassifier(
        **params,
        max_cat_to_onehot=1,  # Helps with categorical features
        #use_label_encoder=False  # Saves memory, deprecated parameters have been removed
    )
    
    print(f"Model will use: {model.get_params().get('nthread', 'default')} threads")
    if has_gpu:
        print(f"GPU is ENABLED - training will be accelerated")
    else:
        print(f"GPU is DISABLED - using CPU only")
    
    # Monitor memory during training
    print("\nTraining started - monitoring memory...")
    
    # Setup memory monitoring thread (optional, for info only)
    import threading
    memory_log = []
    
    def log_memory():
        while training_in_progress:
            mem = psutil.virtual_memory()
            memory_log.append((time.time(), mem.percent))
            time.sleep(10)
    
    training_in_progress = True
    monitor_thread = threading.Thread(target=log_memory, daemon=True)
    monitor_thread.start()
    
    start_time = time.time()
    
    # Fit the model with memory-efficient settings
    model.fit(
        X_train_np, y_train_np,
        eval_set=[(X_train_np, y_train_np), (X_val_np, y_val_np)],
        verbose=50
    )
    
    training_in_progress = False
    training_time = time.time() - start_time
    
    # Memory monitoring results
    if memory_log:
        peak_memory = max(m[1] for m in memory_log)
        print(f"\nPeak memory usage during training: {peak_memory}%")
    
    # ========== MODIFIED: CPU usage monitoring (will be lower with GPU) ==========
    cpu_percent = psutil.cpu_percent(interval=1, percpu=True)
    print(f"\nCPU usage per core: {cpu_percent}")
    print(f"Average CPU usage: {sum(cpu_percent)/len(cpu_percent):.1f}%")
    print(f"Peak threads used: {psutil.Process().num_threads()}")
    print(f"Training time: {training_time/60:.1f} minutes")
    
    # ADDED: GPU usage summary
    if has_gpu:
        try:
            import subprocess
            result = subprocess.run(['nvidia-smi', '--query-gpu=utilization.gpu', '--format=csv,noheader'],
                                  capture_output=True, text=True)
            print(f"GPU utilization during training: {result.stdout.strip()}")
        except:
            pass
    # =======================================================================
    print("="*50)

    # Save the model
    joblib.dump(model, f'{model_dir}/pi0_classifier_model_{br_type}.pkl', compress=3)
    print(f"Model saved to {model_dir}/pi0_classifier_model_{br_type}.pkl")

    print(f"\n{'='*60}")
    print(f"STEP 3: SAVE METRICS")
    print(f"{'='*60}")
    
    # Save metrics
    from sklearn.metrics import roc_auc_score
    
    # Use predict_proba on validation set
    y_pred = model.predict_proba(X_val_np)[:, 1]
    auc_score = roc_auc_score(y_val_np, y_pred)
            
    metrics = {
         'auc': float(auc_score),
         'best_iteration': model.best_iteration,
         'best_score': model.best_score,
         'params': params,
         'n_features': len(training_cols),
         'feature_names': training_cols,
         'training_time_minutes': training_time / 60,
         'peak_memory_percent': peak_memory if memory_log else None,
         'gpu_enabled': has_gpu  # ADDED: track if GPU was used
         }
            
    with open(f'{model_dir}/metrics_{br_type}.json', 'w') as f:
        json.dump(metrics, f, indent=2)
        print(f"Metrics saved to {model_dir}/metrics_{br_type}.json")
        print(f"AUC: {auc_score:.4f}")

    # Get the booster
    booster = model.get_booster()

    # Import the module and patch it
    import ROOT._pythonization._tmva._tree_inference as tree_inference
    tree_inference.get_basescore = patched_get_basescore

    # Save to ROOT
    try:
        ROOT.TMVA.Experimental.SaveXGBoost(
            model,
            "BDT_pi0", 
            f"{model_dir}/bdt_pi0_{br_type}.root", 
            num_inputs=X_train_np.shape[1]
        )
        print(f"✓ Model saved to {model_dir}/bdt_pi0_{br_type}.root")
    except Exception as e:
        print(f"✗ Failed to save ROOT model: {e}")
        # Try alternative save method
        try:
            ROOT.TMVA.Experimental.SaveXGBoost(
                booster,
                "BDT_pi0", 
                f"{model_dir}/bdt_pi0_{br_type}.root", 
                num_inputs=X_train_np.shape[1]
            )
            print(f"✓ Model saved to {model_dir}/bdt_pi0_{br_type}.root (using booster)")
        except Exception as e2:
            print(f"✗ Alternative also failed: {e2}")

    # Memory cleanup
    del X_train_np, X_val_np, y_train_np, y_val_np
    free_memory()
            
    print(f"\n✓ Training completed for {br_type}")