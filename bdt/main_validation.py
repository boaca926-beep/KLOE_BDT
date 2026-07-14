# Validation script
import os
import sys
import joblib
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import xgboost as xgb
from plots import plot_learning_curves, plot_learning_curves_improved, plot_roc, plot_cm, plot_cm_improved, plot_var_score, plot_mass_signal_breakdown, plot_score_breakdown
from metrics import eval_performance, event_performance
from sklearn.metrics import roc_auc_score, roc_curve, auc

import matplotlib
matplotlib.use('Agg')  # Changed from 'TkAgg' to 'Agg' for non-interactive
import matplotlib.pyplot as plt
# plt.show(block=False)  # Commented out
import pandas as pd


from config import (
    DATA_DIR, 
    PLOT_VAL_DIR,
    MODEL_DIR
)

# uv run main_validation.py 2>&1 | tee validation_log.txt

# ========== ADD THIS LINE - CHOOSE DATASET ==========
USE_TEST_SET = False  # Set to False for validation, True for test set
# ====================================================

def load_data(br_nm, input_data_dir):  # ADDED parameters
    """
    Load validation data and models
    """
    
    # ========== MODIFIED: Choose between val and test ==========
    if USE_TEST_SET:
        data_type = 'test'
        print(f"\n📊 Loading TEST set (20% hold-out)...")
    else:
        data_type = 'val'
        print(f"\n📊 Loading VALIDATION set...")
    # ============================================================

    # ADDED: Convert to parquet on first run
    parquet_x = os.path.join(input_data_dir, f'X_{data_type}_{br_nm}.parquet')
    if not os.path.exists(parquet_x):
        print("Converting to parquet (one-time conversion for faster future loads)...")
        for f in [f'X_{data_type}', f'y_{data_type}', f'all_df_{data_type}']:
            pkl_file = os.path.join(input_data_dir, f'{f}_{br_nm}.pkl')
            parquet_file = os.path.join(input_data_dir, f'{f}_{br_nm}.parquet')
            if os.path.exists(pkl_file) and not os.path.exists(parquet_file):
                try:
                    print(f"  Reading {f}_{br_nm}.pkl...")
                    data = joblib.load(pkl_file)
                        
                    if not isinstance(data, pd.DataFrame):
                        if isinstance(data, pd.Series):
                            data = data.to_frame()
                        else:
                            data = pd.DataFrame(data)
                        
                    data.to_parquet(parquet_file, compression='snappy')
                        
                    old_size = os.path.getsize(pkl_file) / 1024**2
                    new_size = os.path.getsize(parquet_file) / 1024**2
                    print(f"  Converted {f}_{br_nm}.pkl: {old_size:.1f}MB → {new_size:.1f}MB")
                except Exception as e:
                    print(f"  Error converting {f}_{br_nm}.pkl: {e}")

    # Load parquet if exists, otherwise load pickle
    if os.path.exists(parquet_x):
        X = pd.read_parquet(os.path.join(input_data_dir, f'X_{data_type}_{br_nm}.parquet'))
        y = pd.read_parquet(os.path.join(input_data_dir, f'y_{data_type}_{br_nm}.parquet'))
        all_df = pd.read_parquet(os.path.join(input_data_dir, f'all_df_{data_type}_{br_nm}.parquet'))
        if isinstance(y, pd.DataFrame):
            y = y.iloc[:, 0]
    else:
        X = joblib.load(os.path.join(input_data_dir, f'X_{data_type}_{br_nm}.pkl'))
        y = joblib.load(os.path.join(input_data_dir, f'y_{data_type}_{br_nm}.pkl'))
        all_df = joblib.load(os.path.join(input_data_dir, f'all_df_{data_type}_{br_nm}.pkl'))
    
    print(f"Loaded {len(X)} events, signal fraction: {y.mean():.4f}")
    
    return X, y, all_df

if __name__ == '__main__':
    print(f"Validation ...")

    input_data_dir = DATA_DIR
    phys_map = joblib.load(os.path.join(input_data_dir, f'phys_map.pkl'))
    
    #print(phys_map)

    ## Load dataset
    phys_ch = ['TCOMB', 'combined']
    br_nm = phys_ch[0]
    info = phys_map.get(br_nm, "")
    #print(info)
    br_title = info['br_title']

    X_val, y_val, all_df = load_data(br_nm, input_data_dir)  # FIXED: Added parameters
    model = joblib.load(os.path.join(MODEL_DIR, f'pi0_classifier_model_{br_nm}.pkl'))
    #print(f"Loading model from: {model}")
    
    plot_dir = PLOT_VAL_DIR
    import shutil
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir, exist_ok=True)

    features = X_val.columns

    ## Evaluate validation set
    eval_performance(model, X_val, y_val)

    ## Feature importance
    importance = model.feature_importances_
    for f, imp in zip(features, importance):
        print(f"    {f}: {imp:.03f}")

    ## 1. Learning curves
    fig_learning = plot_learning_curves_improved(model, rf'Validation')
    #fig_learning = plot_learning_curves(model, rf'Learning Curve (validation)')
    fig_learning.savefig(f'{plot_dir}/learning_curves_{br_nm}.png', dpi=300, bbox_inches='tight')
    plt.close(fig_learning)

    ## 2. Plot confusion matrix
    fig_cm = plot_cm_improved(X_val, y_val, model, rf'Validation')
    #fig_cm = plot_cm(X_val, y_val, model, rf'Confusion Matrix (validation)')
    fig_cm.savefig(f'{plot_dir}/cm_{br_nm}.png', dpi=300, bbox_inches='tight')
    plt.close(fig_cm)
        
    ## 3. Accuracy metrics, event basis
    score_list, var_list, var_str = event_performance(all_df, model)
    #fig_var = plot_var_score(var_list, score_list, var_str, rf'Variable Performance - {br_title}')
    #fig_var.savefig(f'{plot_dir}/var_score_{br_nm}.png', dpi=300, bbox_inches='tight')
    #plt.close(fig_var)

    ## 4. Mass breakdown plot
    y_pred = model.predict(X_val)
    fig_mass = plot_mass_signal_breakdown(X_val, y_val, y_pred, 
                                      phys_ch="TCOMB", 
                                      plot_title="TCOMB: π⁰ Mass – Good Signal vs Lost Signal vs Background")
    fig_mass.savefig(f'{plot_dir}/mass_breakdown_TCOMB.png', dpi=300, bbox_inches='tight')
    plt.close(fig_mass)

    # 5. Score breakdown
    y_score = model.predict_proba(X_val)[:, 1]
    fig_score = plot_score_breakdown(y_val, y_pred, y_score,
                                    phys_ch=br_nm,
                                    plot_title=f"{br_title}: BDT Score – Signal vs Background")
    fig_score.savefig(f"{plot_dir}/score_breakdown_{br_nm}.png", dpi=300, bbox_inches='tight')
    plt.close(fig_score)

    ## 6. ROC plot
    #fig_roc = plot_roc(score_list, rf'ROC Curve - $\pi^{0}$ Classifier (validation, {br_title})')
    #fig_roc.savefig(f'{plot_dir}/roc_curv_{br_nm}.png', dpi=300, bbox_inches='tight')
    #plt.close(fig_roc)

    y_score = model.predict_proba(X_val)[:, 1]
    #fig_roc = plot_roc(y_val, y_score, rf'ROC Curve - π⁰ Classifier (validation)')
    #fig_roc.savefig(f'{plot_dir}/roc_curv_{br_nm}.png', dpi=300, bbox_inches='tight')
    #plt.close(fig_roc)

