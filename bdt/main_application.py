import pandas as pd
import numpy as np
import os
import sys
import joblib
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from plots import plot_var_score, plot_roc, plot_nm  # Fixed import
from metrics import event_performance
from sklearn.metrics import confusion_matrix
import seaborn as sns

from config import DATA_DIR, PLOT_APP_DIR, MODEL_DIR

def event_wise_prediction(all_df_test, X_test, y_test_pair, model, threshold=0.5):
    """
    Convert pair-wise predictions to event-wise decisions
    """
    
    # Get pair predictions from model
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    y_pred_pair = (y_pred_proba >= threshold).astype(int)
    
    # ============ CRITICAL FIX: Verify data alignment ============
    n_photons = len(all_df_test)
    n_pairs = len(y_pred_proba)
    n_pairs_per_photon = 3
    
    # Check if we have the expected number of pairs
    expected_pairs = n_photons * n_pairs_per_photon
    if n_pairs != expected_pairs:
        print(f"⚠️ WARNING: Expected {expected_pairs} pairs, got {n_pairs}")
        print(f"   Using available pairs only")
        n_photons = min(n_photons, n_pairs // n_pairs_per_photon)
    
    print(f"\nData validation:")
    print(f"  Photons: {n_photons}")
    print(f"  Pairs: {n_pairs}")
    print(f"  Pairs per photon: {n_pairs_per_photon}")
    
    # ============ PHOTON-LEVEL AGGREGATION ============
    photon_data = []
    
    for photon_idx in range(n_photons):
        start_idx = photon_idx * n_pairs_per_photon
        end_idx = start_idx + n_pairs_per_photon
        
        # Safety check
        if end_idx > n_pairs:
            print(f"  Truncating at photon {photon_idx}")
            break
        
        # Get photon info
        event_id = all_df_test.iloc[photon_idx]['event']
        is_signal_photon = all_df_test.iloc[photon_idx].get('is_signal', 0)
        
        # Get predictions for this photon's 3 pairs
        photon_probs = y_pred_proba[start_idx:end_idx]
        photon_preds = y_pred_pair[start_idx:end_idx]
        
        # Check if this photon is part of a true π⁰ pair
        has_true_pi0 = False
        true_pair_count = 0
        for pair_offset in range(n_pairs_per_photon):
            pair_idx = start_idx + pair_offset
            if pair_idx < len(y_test_pair) and y_test_pair.iloc[pair_idx] == 1:
                has_true_pi0 = True
                true_pair_count += 1
        
        photon_data.append({
            'photon_idx': photon_idx,
            'event_id': event_id,
            'is_signal_photon': is_signal_photon,
            'has_true_pi0': has_true_pi0,
            'true_pair_count': true_pair_count,
            'max_proba': photon_probs.max(),
            'mean_proba': photon_probs.mean(),
            'n_pred_pi0_pairs': photon_preds.sum(),
            'has_pred_pi0_pair': int(photon_preds.sum() > 0),
        })
    
    photon_df = pd.DataFrame(photon_data)
    
    # ============ EVENT-LEVEL AGGREGATION ============
    event_results = []
    
    for event_id, group in photon_df.groupby('event_id'):
        n_photons_in_event = len(group)
        
        # TRUE LABEL: Based on pair-level truth
        true_signal = int(group['has_true_pi0'].sum() > 0)
        
        # Track how many true π⁰ pairs in this event
        total_true_pairs = group['true_pair_count'].sum()
        
        # PREDICTIONS:
        pred_any = int(group['has_pred_pi0_pair'].sum() > 0)
        pred_min2 = int(group['has_pred_pi0_pair'].sum() >= 2)
        pred_max = int(group['max_proba'].max() >= threshold)
        pred_mean = int(group['mean_proba'].mean() >= threshold)
        
        event_results.append({
            'event_id': event_id,
            'true_signal': true_signal,
            'total_true_pairs': total_true_pairs,
            'n_photons': n_photons_in_event,
            'n_signal_photons': group['is_signal_photon'].sum(),
            'n_photons_with_true_pi0': group['has_true_pi0'].sum(),
            'n_photons_with_pred_pi0': group['has_pred_pi0_pair'].sum(),
            'max_proba': group['max_proba'].max(),
            'mean_proba': group['mean_proba'].mean(),
            'pred_any': pred_any,
            'pred_min2': pred_min2,
            'pred_max': pred_max,
            'pred_mean': pred_mean
        })
    
    event_df = pd.DataFrame(event_results)
    
    # ============ EVALUATION ============
    print("\n" + "="*60)
    print("EVENT-LEVEL ANALYSIS")
    print("="*60)
    
    print(f"\nTrue event distribution (based on π⁰ pairs):")
    true_signal_count = event_df['true_signal'].sum()
    print(f"  Signal events: {true_signal_count}")
    print(f"  Background events: {len(event_df) - true_signal_count}")
    
    # Evaluate strategies
    strategies = [
        ('any', 'pred_any'),
        ('min2', 'pred_min2'),
        ('max', 'pred_max'),
        ('mean', 'pred_mean')
    ]
    
    print(f"\nStrategy Performance (threshold={threshold}):")
    best_f1 = 0
    best_strategy = 'any'
    
    for strategy_name, col in strategies:
        tp = ((event_df['true_signal'] == 1) & (event_df[col] == 1)).sum()
        fp = ((event_df['true_signal'] == 0) & (event_df[col] == 1)).sum()
        tn = ((event_df['true_signal'] == 0) & (event_df[col] == 0)).sum()
        fn = ((event_df['true_signal'] == 1) & (event_df[col] == 0)).sum()
        
        if (tp + fp + tn + fn) > 0:
            accuracy = (tp + tn) / len(event_df)
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
            
            print(f"\n  Strategy '{strategy_name}':")
            print(f"    TP: {tp:6d} | FP: {fp:6d} | TN: {tn:6d} | FN: {fn:6d}")
            print(f"    Precision: {precision:.4f}, Recall: {recall:.4f}, F1: {f1:.4f}")
            
            if f1 > best_f1:
                best_f1 = f1
                best_strategy = strategy_name
    
    print(f"\n✓ Best strategy at threshold={threshold}: '{best_strategy}' (F1 = {best_f1:.4f})")
    
    # Add the best strategy as default for this threshold
    event_df['pred_signal'] = event_df[f'pred_{best_strategy}']
    
    return event_df, best_strategy, best_f1


def plot_event_confusion_matrix(event_results, data_type, plot_dir):
    """Plot confusion matrix for event-wise classification"""
    
    cm = confusion_matrix(event_results['true_signal'], event_results['pred_signal'])
    cm_percent = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot counts
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0],
                xticklabels=['Background', 'Signal'],
                yticklabels=['Background', 'Signal'])
    axes[0].set_xlabel('Predicted')
    axes[0].set_ylabel('True')
    axes[0].set_title(f'Confusion Matrix (Counts) - {data_type}')
    
    # Plot percentages
    sns.heatmap(cm_percent, annot=True, fmt='.1f', cmap='Blues', ax=axes[1],
                xticklabels=['Background', 'Signal'],
                yticklabels=['Background', 'Signal'])
    axes[1].set_xlabel('Predicted')
    axes[1].set_ylabel('True')
    axes[1].set_title(f'Confusion Matrix (Percentages) - {data_type}')
    
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/event_cm_{data_type}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Calculate and print metrics
    tn, fp, fn, tp = cm.ravel()
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\nEvent-wise Classification Metrics ({data_type}):")
    print(f"  True Positives:  {tp:5d}")
    print(f"  False Positives: {fp:5d}")
    print(f"  True Negatives:  {tn:5d}")
    print(f"  False Negatives: {fn:5d}")
    print(f"  Accuracy:  {accuracy:.4f}")
    print(f"  Precision: {precision:.4f}")
    print(f"  Recall:    {recall:.4f}")
    print(f"  F1 Score:  {f1:.4f}")
    
    return fig


if __name__ == '__main__':
    
    print(f"Application on test dataset...")
    
    input_data_dir = DATA_DIR
    input_model_dir = MODEL_DIR
    
    category_type = 'TCOMB'
    
    # Create output folder
    plot_dir = PLOT_APP_DIR
    
    import shutil
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir, exist_ok=True)
    
    # Load physics map
    phys_map = joblib.load(os.path.join(input_data_dir, f'phys_map.pkl'))
    
    for data_type, info in phys_map.items():
        br_title = info['br_title']
        category = info['category']
        
        print(f"\n{'='*60}")
        print(f"Processing: {data_type}")
        print(f"{'='*60}")
        
        if data_type == category_type:
            
            # Load data
            all_df = joblib.load(os.path.join(input_data_dir, f'all_df_{data_type}.pkl'))
            all_df_test = joblib.load(os.path.join(input_data_dir, f'all_df_test_{data_type}.pkl'))
            X_test = joblib.load(os.path.join(input_data_dir, f'X_test_{data_type}.pkl'))
            y_test = joblib.load(os.path.join(input_data_dir, f'y_test_{data_type}.pkl'))
            model = joblib.load(os.path.join(input_model_dir, f'pi0_classifier_model_{data_type}.pkl'))
            
            print(f"Loaded: {len(all_df_test)} photons, {len(X_test)} pairs")
            
            # Get event-wise prediction
            #event_results, best_strategy, best_f1 = event_wise_prediction(
            #    all_df_test, X_test, y_test, model, threshold=0.5
            #)

            # Analyze threshold impact
            thresholds = np.arange(0.05, 1.0, 0.05)
            best_f1 = 0.0
            best_thr = 0.5
            best_strat = 'any'
            best_event_df = None

            for thr in thresholds:
                event_df, strat, f1 = event_wise_prediction(
                all_df_test, X_test, y_test, model, threshold=thr
                )
                if f1 > best_f1:
                    best_f1 = f1
                    best_thr = thr
                    best_strat = strat
                    best_event_df = event_df

            print(f"\n{'='*60}")
            print(f"OPTIMAL CONFIGURATION")
            print(f"{'='*60}")
            print(f"Threshold  : {best_thr:.2f}")
            print(f"Strategy   : '{best_strat}'")
            print(f"F1 score   : {best_f1:.4f}")

            # Use best_event_df for plotting and saving
            event_results = best_event_df

            
            # Plot event confusion matrix
            plot_event_confusion_matrix(event_results, data_type, plot_dir)
            
            # Save results
            event_results.to_csv(f'{plot_dir}/event_results_{data_type}.csv', index=False)
            
            # Plot pair-level metrics
            fig_cm = plot_nm(X_test, y_test, model, br_title)
            fig_cm.savefig(f'{plot_dir}/cm_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_cm)
            
            # Plot ROC and mass-score
            score_list, var_list, var_str = event_performance(all_df, model)
            fig_roc = plot_roc(score_list, rf'ROC Curve - $\pi^{0}$ Classifier (test, {br_title})')
            fig_roc.savefig(f'{plot_dir}/roc_curv_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_roc)
            
            print(f"\n✓ All results saved to {plot_dir}")
        
        else:
            print(f"Skipping {data_type} - not target category")