import pandas as pd
import numpy as np
import os
import sys
import joblib
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from plots import plot_roc_improved, plot_event_cm_improved, plot_cm_improved, plot_f1_vs_threshold
from sklearn.metrics import roc_curve, auc
from config import DATA_DIR, PLOT_APP_DIR, MODEL_DIR

# uv run main_application.py 2>&1 | tee application_log.txt


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def aggregate_event_scores(all_df_test, y_pred_proba, y_test_pair, pairs_per_photon=3):
    """
    Convert pair-level probabilities to event-level continuous scores.
    
    Returns a DataFrame with columns:
        event_id, true_signal, max_score, mean_score, min2_score
    
    - max_score: maximum pair probability in the event (for 'any' strategy)
    - mean_score: average of all pair probabilities in the event (for 'mean' strategy)
    - min2_score: second highest per‑photon maximum (for 'min2' strategy)
    """
    n_photons = len(all_df_test)
    n_pairs = len(y_pred_proba)
    
    # Validate alignment
    expected_pairs = n_photons * pairs_per_photon
    assert n_pairs == expected_pairs, \
        f"Mismatch: {n_photons} photons × {pairs_per_photon} ≠ {n_pairs} pairs"
    
    y_test_array = y_test_pair.values if hasattr(y_test_pair, 'values') else y_test_pair
    
    photon_records = []
    for i in range(n_photons):
        start = i * pairs_per_photon
        end = start + pairs_per_photon
        probs = y_pred_proba[start:end]
        labels = y_test_array[start:end]
        
        photon_records.append({
            'event_id': all_df_test.iloc[i]['event'],
            'photon_max': probs.max(),          # max within this photon
            'photon_mean': probs.mean(),        # mean within this photon
            'has_true_pi0': labels.sum() > 0    # photon contains a true π⁰ pair
        })
    
    photon_df = pd.DataFrame(photon_records)
    
    # Define helper to compute second largest value
    def second_largest(arr):
        sorted_arr = np.sort(arr)
        if len(sorted_arr) >= 2:
            return sorted_arr[-2]
        else:
            # Events with <2 photons cannot be signal under min2
            return 0.0
    
    # Aggregate at event level
    event_df = photon_df.groupby('event_id').agg({
        'photon_max': 'max',                    # overall maximum -> any strategy
        'photon_mean': 'mean',                  # overall mean -> mean strategy
        'has_true_pi0': 'any'                   # event truth
    }).rename(columns={
        'photon_max': 'max_score',
        'photon_mean': 'mean_score',
        'has_true_pi0': 'true_signal'
    }).reset_index()
    
    # Compute min2_score: second highest of per‑photon maxima
    event_df['min2_score'] = photon_df.groupby('event_id')['photon_max'].agg(second_largest).values
    
    # Convert boolean truth to int
    event_df['true_signal'] = event_df['true_signal'].astype(int)
    
    return event_df


def event_wise_prediction_fast(all_df_test, y_pred_proba, y_test_pair, threshold=0.5):
    """
    Convert pair-wise predictions to event-wise decisions.
    Uses precomputed probabilities to avoid repeated model inference.
    Returns event_df, best_strategy, best_f1.
    """
    # Get pair predictions from precomputed probabilities
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
    f1_any = 0.0
    
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
            
            # fill f1_any
            if (strategy_name == 'any'):
                f1_any = f1

            print(f"\n  Strategy '{strategy_name}':")
            print(f"    Accuracy '{accuracy:.4f}'")
            print(f"    TP: {tp:6d} | FP: {fp:6d} | TN: {tn:6d} | FN: {fn:6d}")
            print(f"    Precision: {precision:.4f}, Recall: {recall:.4f}, F1: {f1:.4f}")
            
            if f1 > best_f1:
                best_f1 = f1
                best_strategy = strategy_name
    
    print(f"\n✓ Best strategy at threshold={threshold}: '{best_strategy}' (F1 = {best_f1:.4f})")
    #print(f"\n f1 'any' list {f1_any}")

    # Add the best strategy as default for this threshold
    event_df['pred_signal'] = event_df[f'pred_{best_strategy}']
    
    return event_df, best_strategy, best_f1, f1_any


# =============================================================================
# MAIN EXECUTION
# =============================================================================

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
            
            # Validate alignment
            assert len(all_df_test) * 3 == len(X_test), \
                f"Mismatch: {len(all_df_test)} photons × 3 ≠ {len(X_test)} pairs"
            
            # ========== OPTIMIZATION: Compute probabilities once ==========
            print("\nComputing pair-level probabilities (once)...")
            y_pred_proba = model.predict_proba(X_test)[:, 1]
            print(f"Done. Probabilities shape: {y_pred_proba.shape}")
            
            # =====================================================================
            # EVENT-LEVEL ROC CURVES – using event aggregation
            # =====================================================================
            print("\n" + "="*60)
            print("AGGREGATING TO EVENT-LEVEL SCORES FOR ROC")
            print("="*60)
            
            event_scores = aggregate_event_scores(all_df_test, y_pred_proba, y_test)
            print(f"✅ Aggregated {len(event_scores)} events.")
            
            # ---- Plot and save event-wise ROC curves ----
            # 1. Max probability (any strategy)
            fig_roc_max = plot_roc_improved(
                event_scores['true_signal'], event_scores['max_score'],
                plot_title=f'ROC Curve - π⁰ Classifier (max probability)',
                threshold=0.35
            )
            fig_roc_max.savefig(f'{plot_dir}/event_roc_max_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_roc_max)
            
            # 2. Mean probability (mean strategy)
            fig_roc_mean = plot_roc_improved(
                event_scores['true_signal'], event_scores['mean_score'],
                plot_title=f'Event-wise ROC (mean probability) – {data_type}'
            )
            fig_roc_mean.savefig(f'{plot_dir}/event_roc_mean_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_roc_mean)
            
            # 3. Second-highest photon-max (min2 strategy)
            fig_roc_second = plot_roc_improved(
                event_scores['true_signal'], event_scores['min2_score'],
                plot_title=f'Event-wise ROC (2nd max photon probability) – {data_type}'
            )
            fig_roc_second.savefig(f'{plot_dir}/event_roc_min2_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_roc_second)
            
            # ---- Optional: overlay all three strategies ----
            fig_overlay, ax = plt.subplots(figsize=(10, 8))
            for scores, label, color in zip(
                [event_scores['max_score'], event_scores['mean_score'], event_scores['min2_score']],
                ['Max (any)', 'Mean', '2nd Max (min2)'],
                ['#1f77b4', '#ff7f0e', '#2ca02c']
            ):
                fpr, tpr, _ = roc_curve(event_scores['true_signal'], scores)
                auc_val = auc(fpr, tpr)
                ax.plot(fpr, tpr, lw=2, color=color, label=f'{label} (AUC={auc_val:.3f})')
            ax.plot([0, 1], [0, 1], 'k--', lw=1, label='Random')
            ax.set_xlabel('False Positive Rate', fontsize=14)
            ax.set_ylabel('True Positive Rate', fontsize=14)
            ax.set_title(f'Event-wise ROC Comparison – {data_type}', fontsize=16)
            ax.legend(loc='lower right', fontsize=12)
            ax.grid(True, alpha=0.3)
            fig_overlay.savefig(f'{plot_dir}/event_roc_comparison_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_overlay)
            
            print(f"✅ Event-wise ROC curves saved to: {plot_dir}")
            
            # =====================================================================
            # THRESHOLD OPTIMISATION (discrete strategies)
            # =====================================================================
            thresholds = np.arange(0.05, 1, 0.05)
            best_f1 = 0.0
            best_thr = 0.5
            best_strat = 'any'
            best_event_df = None
            
            all_f1_any = []
            thr_list = []

            for thr in thresholds:
                event_df, strat, f1, f1_any = event_wise_prediction_fast(
                    all_df_test, y_pred_proba, y_test, threshold=thr
                )
                all_f1_any.append(f1_any)
                thr_list.append(thr)
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
            
            # Use best_event_df for saving and confusion matrix
            event_results = best_event_df
            
            # Save results
            event_results.to_csv(f'{plot_dir}/event_results_{data_type}.csv', index=False)
            
            # ---- Plot event confusion matrix ----
            fig_event_cm = plot_event_cm_improved(event_results, data_type)
            fig_event_cm.savefig(f'{plot_dir}/event_cm_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_event_cm)
            
            # ---- (Optional) Pair-level confusion matrix ----
            fig_cm = plot_cm_improved(X_test, y_test, model, br_title)
            fig_cm.savefig(f'{plot_dir}/cm_pair_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_cm)
            
            # ---- (Optional) Pair-level ROC ----
            fig_pair_roc = plot_roc_improved(y_test, y_pred_proba, 
                                             plot_title=f'Pair-level ROC – {data_type}')
            fig_pair_roc.savefig(f'{plot_dir}/event_roc_pair_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_pair_roc)
            
            # F1 plot 'any'
            clean_f1 = [float(f"{x:.4f}") for x in all_f1_any]
            print(f"all_f1_any = {all_f1_any}")

            fig_f1_any = plot_f1_vs_threshold(clean_f1, thr_list)   # now works with default opt_threshold
            fig_f1_any.savefig(f'{plot_dir}/f1_any.png', dpi=300, bbox_inches='tight')
            plt.close(fig_f1_any)

            print(f"\n✓ All results saved to {plot_dir}")
        
        else:
            print(f"Skipping {data_type} - not target category")