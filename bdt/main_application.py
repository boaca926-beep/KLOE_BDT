import pandas as pd
import numpy as np
import os
import sys
import joblib
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from plots import plot_roc_improved, plot_event_cm_improved, plot_cm_improved, plot_f1_vs_threshold, plot_event_score_breakdown
from sklearn.metrics import roc_curve, auc
from config import DATA_DIR, PLOT_APP_DIR, MODEL_DIR

# uv run main_application.py 2>&1 | tee application_log.txt


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def aggregate_event_scores(all_df_test, y_pred_proba, y_test_pair, pairs_per_event=3):
    """
    Convert pair-level probabilities to event-level continuous scores.
    Assumes all_df_test has one row per event, and y_pred_proba has 3 probabilities per event.
    """
    n_events = len(all_df_test)
    n_pairs = len(y_pred_proba)
    
    # Validate: 3 pairs per event
    assert n_pairs == n_events * pairs_per_event, \
        f"Mismatch: {n_events} events × {pairs_per_event} ≠ {n_pairs} pairs"
    
    # Reshape probabilities to (n_events, 3)
    probs_matrix = y_pred_proba.reshape(n_events, pairs_per_event)
    y_test_matrix = y_test_pair.values.reshape(n_events, pairs_per_event)
    
    # Compute scores
    max_score = probs_matrix.max(axis=1)
    mean_score = probs_matrix.mean(axis=1)
    
    # Second largest (for min2)
    sorted_probs = np.sort(probs_matrix, axis=1)
    min2_score = sorted_probs[:, 1]   # 2nd column (0-indexed) -> second largest
    
    # Event truth: any of the 3 pairs is signal
    true_signal = (y_test_matrix.sum(axis=1) > 0).astype(int)
    
    event_df = pd.DataFrame({
        'event_id': all_df_test['event'].values,
        'true_signal': true_signal,
        'max_score': max_score,
        'mean_score': mean_score,
        'min2_score': min2_score
    })
    
    return event_df


def event_wise_prediction_fast(all_df_test, y_pred_proba, y_test_pair, threshold=0.5):
    """
    Convert pair-wise predictions to event-wise decisions.
    all_df_test: one row per event.
    y_pred_proba: 3 probabilities per event (concatenated).
    y_test_pair: 3 labels per event (concatenated).
    Returns: event_df, best_strategy, best_f1, f1_scores_dict
    """
    n_events = len(all_df_test)
    n_pairs = len(y_pred_proba)
    pairs_per_event = 3
    
    # Safety check
    if n_pairs != n_events * pairs_per_event:
        print(f"⚠️ WARNING: Expected {n_events * pairs_per_event} pairs, got {n_pairs}")
        n_events = n_pairs // pairs_per_event
        print(f"   Truncating to {n_events} events")
    
    # Reshape
    probs = y_pred_proba[:n_events * pairs_per_event].reshape(n_events, pairs_per_event)
    labels = y_test_pair.values[:n_events * pairs_per_event].reshape(n_events, pairs_per_event)
    
    event_results = []
    for i in range(n_events):
        event_id = all_df_test.iloc[i]['event']
        p = probs[i]
        lbl = labels[i]
        
        true_signal = int(lbl.sum() > 0)
        pred_any = int(p.max() >= threshold)
        pred_min2 = int(np.sort(p)[1] >= threshold)   # second largest
        pred_mean = int(p.mean() >= threshold)
        # pred_max is the same as pred_any for this setup
        pred_max = pred_any
        
        event_results.append({
            'event_id': event_id,
            'true_signal': true_signal,
            'max_proba': p.max(),
            'mean_proba': p.mean(),
            'min2_proba': np.sort(p)[1],
            'pred_any': pred_any,
            'pred_min2': pred_min2,
            'pred_max': pred_max,
            'pred_mean': pred_mean
        })
    
    event_df = pd.DataFrame(event_results)
    
    # Evaluate strategies and pick best F1
    strategies = [('any', 'pred_any'), ('min2', 'pred_min2'), 
                  ('max', 'pred_max'), ('mean', 'pred_mean')]
    
    best_f1 = 0
    best_strategy = 'any'
    f1_scores = {}   # store per-strategy F1
    
    for name, col in strategies:
        tp = ((event_df['true_signal'] == 1) & (event_df[col] == 1)).sum()
        fp = ((event_df['true_signal'] == 0) & (event_df[col] == 1)).sum()
        fn = ((event_df['true_signal'] == 1) & (event_df[col] == 0)).sum()
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        f1_scores[name] = f1
        
        if f1 > best_f1:
            best_f1 = f1
            best_strategy = name
    
    event_df['pred_signal'] = event_df[f'pred_{best_strategy}']
    
    return event_df, best_strategy, best_f1, f1_scores


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
                plot_title=f'Event-wise ROC (max probability) – {data_type}'
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
                #[event_scores['max_score'], event_scores['mean_score'], event_scores['min2_score']],
                #['Max (any)', 'Mean', '2nd Max (min2)'],
                #['#1f77b4', '#ff7f0e', '#2ca02c']
                [event_scores['max_score'], event_scores['mean_score']],
                ['Max (any)', 'Mean'],
                ['#1f77b4', '#ff7f0e']
            ):
                fpr, tpr, _ = roc_curve(event_scores['true_signal'], scores)
                auc_val = auc(fpr, tpr)
                print(f"auc_val: {label}={auc_val:.3}")
                ax.plot(fpr, tpr, lw=2, color=color, label=f'{label} (AUC={auc_val:.3f})')
            ax.plot([0, 1], [0, 1], 'k--', lw=1, label='Random')
            ax.set_xlabel('False Positive Rate', fontsize=20)
            ax.set_ylabel('True Positive Rate', fontsize=20)
            ax.tick_params(axis='both', labelsize=18)
            title=f'Event-wise ROC Comparison – {data_type}'
            ax.set_title('', fontsize=16)
            ax.legend(loc='lower right', fontsize=20)
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
            all_f1_mean = []
            all_f1_min2 = []   # optional
            thr_list = []

            for thr in thresholds:
                event_df, strat, f1, f1_scores = event_wise_prediction_fast(
                    all_df_test, y_pred_proba, y_test, threshold=thr
                )
                all_f1_any.append(f1_scores['any'])
                all_f1_mean.append(f1_scores['mean'])
                all_f1_min2.append(f1_scores['min2'])   # optional
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
            fig_event_cm, accuracy, precision, recall, f1 = plot_event_cm_improved(event_results, data_type)
            fig_event_cm.savefig(f'{plot_dir}/event_cm_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_event_cm)

            # Print accuracy, precision, recall, f1 (event-wise)
            print(f"\n{'='*60}")
            print(f"Event-wise Metric")
            print(f"{'='*60}")
            print(f"Accuracy  : {accuracy:.4f}")
            print(f"Precision : {precision:.4f}")
            print(f"Recall    : {recall:.4f}")
            print(f"F1 score  : {f1:.4f}")
            
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
            fig_f1_any = plot_f1_vs_threshold(
                all_f1_any, 
                thr_list, 
                strategy_name='any'
            )
            fig_f1_any.savefig(f'{plot_dir}/f1_any.png', dpi=300, bbox_inches='tight')
            plt.close(fig_f1_any)
            print(f"✅ Saved F1 vs threshold for 'any' strategy")
            
            # F1 plot 'mean'
            fig_f1_mean = plot_f1_vs_threshold(
                all_f1_mean,  
                thr_list, 
                strategy_name='mean'
            )
            fig_f1_mean.savefig(f'{plot_dir}/f1_mean.png', dpi=300, bbox_inches='tight')
            plt.close(fig_f1_mean)
            print(f"✅ Saved F1 vs threshold for 'mean' strategy")

            # ---- Event-wise score separation plot ----
            fig_event_score = plot_event_score_breakdown(
                event_results, 
                score_col='mean_proba', # or 'max_proba', 'mean_proba' / 'min2_proba'
                threshold=best_thr,
                phys_ch=data_type
            )
            fig_event_score.savefig(f'{plot_dir}/event_score_breakdown_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_event_score)
            print(f"\n✓ All results saved to {plot_dir}")
        
        else:
            print(f"Skipping {data_type} - not target category")