import pandas as pd
import numpy as np
import os
import sys
import joblib
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from plots import plot_var_score, plot_roc, plot_nm
from metrics import event_performance
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

from config import DATA_DIR, PLOT_APP_DIR, MODEL_DIR

def event_wise_prediction(all_df_test, X_test, y_test_pair, model, threshold=0.5):
    """
    Convert pair-wise predictions to event-wise decisions
    
    FIX: Properly handles the mapping between photons and their pairs
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
        # FIX: Check the actual pair labels for this photon
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
    # FIX: Use the correct truth definition based on your data
    # Since we have events with 0, 1, or 2 signal photons, we need to define
    # what constitutes a signal event based on the pair labels
    
    event_results = []
    
    for event_id, group in photon_df.groupby('event_id'):
        n_photons_in_event = len(group)
        
        # TRUE LABEL: Based on pair-level truth
        # An event is signal if it contains at least one true π⁰ pair
        true_signal = int(group['has_true_pi0'].sum() > 0)
        
        # Track how many true π⁰ pairs in this event
        total_true_pairs = group['true_pair_count'].sum()
        
        # PREDICTIONS:
        # Strategy 1: Event is signal if ANY photon has a predicted π⁰ pair
        pred_any = int(group['has_pred_pi0_pair'].sum() > 0)
        
        # Strategy 2: Event is signal if at least 2 photons have predicted π⁰ pairs
        # (This is more physically motivated for π⁰ → γγ)
        pred_min2 = int(group['has_pred_pi0_pair'].sum() >= 2)
        
        # Strategy 3: Use max probability
        pred_max = int(group['max_proba'].max() >= threshold)
        
        # Strategy 4: Use mean probability
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
    print(f"  Average true π⁰ pairs per signal event: {event_df[event_df['true_signal']==1]['total_true_pairs'].mean():.2f}")
    
    # Check correlation with photon-level signal
    print(f"\nCorrelation with photon-level signal:")
    events_with_2_signal = (event_df['n_signal_photons'] == 2).sum()
    events_with_true_pairs_and_2_signal = ((event_df['n_signal_photons'] == 2) & (event_df['true_signal'] == 1)).sum()
    print(f"  Events with exactly 2 signal photons: {events_with_2_signal}")
    print(f"  ...that are signal by pair definition: {events_with_true_pairs_and_2_signal}")
    
    # Evaluate strategies
    strategies = [
        ('any', 'pred_any'),
        ('min2', 'pred_min2'),
        ('max', 'pred_max'),
        ('mean', 'pred_mean')
    ]
    
    print(f"\nStrategy Performance:")
    best_f1 = 0
    best_strategy = 'any'
    best_threshold_for_strategy = threshold
    
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
                best_threshold_for_strategy = threshold
    
    print(f"\n✓ Best strategy at threshold={threshold}: '{best_strategy}' (F1 = {best_f1:.4f})")
    
    # Add the best strategy as default for this threshold
    event_df['pred_signal'] = event_df[f'pred_{best_strategy}']
    
    return event_df

def plot_event_confusion_matrix(event_results, data_type, plot_dir):
    """
    Plot confusion matrix for event-wise classification
    """
    
    # Use the correct column name 'true_signal' instead of 'true_signal_pair_based'
    # Check if the column exists, if not, use the correct one
    true_col = 'true_signal' if 'true_signal' in event_results.columns else 'true_signal_pair_based'
    pred_col = 'pred_signal' if 'pred_signal' in event_results.columns else 'pred_signal'
    
    # Calculate confusion matrix
    cm = confusion_matrix(event_results[true_col], 
                          event_results[pred_col])
    
    # Calculate percentages
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

def analyze_threshold_impact(event_results, data_type, plot_dir):
    """
    Analyze how different thresholds affect event-wise classification
    """
    thresholds = np.arange(0.05, 1.0, 0.05)
    results = []
    
    for threshold in thresholds:
        row_results = {'threshold': threshold}
        
        # Recalculate predictions with new threshold
        # Note: 'any' cannot be properly recalculated here without per-photon predictions
        # So we use max_proba as a proxy (same as 'max')
        pred_any = (event_results['max_proba'] >= threshold).astype(int)
        pred_max = (event_results['max_proba'] >= threshold).astype(int)
        pred_mean = (event_results['mean_proba'] >= threshold).astype(int)
        pred_min2 = (event_results['n_photons_with_pred_pi0'] >= 2).astype(int)
        
        for strategy, pred in [('any', pred_any), ('min2', pred_min2), 
                                ('max', pred_max), ('mean', pred_mean)]:
            tp = ((event_results['true_signal'] == 1) & (pred == 1)).sum()
            fp = ((event_results['true_signal'] == 0) & (pred == 1)).sum()
            fn = ((event_results['true_signal'] == 1) & (pred == 0)).sum()
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            row_results[f'{strategy}_precision'] = precision
            row_results[f'{strategy}_recall'] = recall
            row_results[f'{strategy}_f1'] = f1
        
        results.append(row_results)
    
    results_df = pd.DataFrame(results)
    
    # Find best strategy and threshold
    best_overall_f1 = 0
    best_strategy = 'any'
    best_threshold = 0.5
    
    strategy_type = ['any', 'min2', 'max', 'mean']
    for strategy in strategy_type:
        f1_col = f'{strategy}_f1'
        if f1_col in results_df.columns:
            max_idx = results_df[f1_col].idxmax()
            best_f1 = results_df.loc[max_idx, f1_col]
            if best_f1 > best_overall_f1:
                best_overall_f1 = best_f1
                best_strategy = strategy
                best_threshold = results_df.loc[max_idx, 'threshold']
    
    # Plot - only 3 subplots (no empty 4th)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    for strategy in strategy_type:
        prec_col = f'{strategy}_precision'
        recall_col = f'{strategy}_recall'
        f1_col = f'{strategy}_f1'
        
        if prec_col in results_df.columns:
            axes[0, 0].plot(results_df['threshold'], results_df[prec_col], linewidth=2, label=strategy)
        if recall_col in results_df.columns:
            axes[0, 1].plot(results_df['threshold'], results_df[recall_col], linewidth=2, label=strategy)
        if f1_col in results_df.columns:
            axes[1, 0].plot(results_df['threshold'], results_df[f1_col], linewidth=2, label=strategy)
    
    axes[0, 0].set_ylabel('Precision')
    axes[0, 1].set_ylabel('Recall')
    axes[1, 0].set_ylabel('F1 Score')
    axes[1, 1].axis('off')
    
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set_xlabel('Threshold')
        ax.grid(True, alpha=0.3)
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend()
    
    plt.suptitle(f'Threshold Impact Analysis - {data_type}')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/threshold_analysis_{data_type}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"\nBest overall: Strategy '{best_strategy}' at threshold {best_threshold:.3f} (F1 = {best_overall_f1:.4f})")
    
    return results_df, best_strategy, best_threshold

if __name__ == '__main__':
    
    print(f"Application on test dataset...")
    
    input_data_dir = DATA_DIR
    input_model_dir = MODEL_DIR
    
    category_type = 'TCOMB'  # Use your category
    
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
            
            # Load data (FIX: Pass directories)
            all_df = joblib.load(os.path.join(input_data_dir, f'all_df_{data_type}.pkl'))
            all_df_test = joblib.load(os.path.join(input_data_dir, f'all_df_test_{data_type}.pkl'))
            X_test = joblib.load(os.path.join(input_data_dir, f'X_test_{data_type}.pkl'))
            y_test = joblib.load(os.path.join(input_data_dir, f'y_test_{data_type}.pkl'))
            model = joblib.load(os.path.join(input_model_dir, f'pi0_classifier_model_{data_type}.pkl'))
            
            print(f"Loaded: {len(all_df_test)} photons, {len(X_test)} pairs")
            
            # Get event-wise prediction
            event_results = event_wise_prediction(
                all_df_test, X_test, y_test, model, threshold=0.5
            )
            
            # Plot event confusion matrix
            plot_event_confusion_matrix(event_results, data_type, plot_dir)
            
            # Analyze threshold impact
            threshold_results, best_strategy, best_threshold = analyze_threshold_impact(
                event_results, data_type, plot_dir
            )
            
            # Re-run with best threshold if needed
            if best_threshold != 0.5:
                print(f"\nRe-running with optimal threshold: {best_threshold:.3f}")
                event_results_optimal = event_wise_prediction(
                    all_df_test, X_test, y_test, model, threshold=best_threshold
                )
                plot_event_confusion_matrix(event_results_optimal, f"{data_type}_opt", plot_dir)
                event_results_optimal.to_csv(f'{plot_dir}/event_results_{data_type}_opt.csv', index=False)
            
            # Save results
            event_results.to_csv(f'{plot_dir}/event_results_{data_type}.csv', index=False)
            
            # Plot pair-level metrics
            fig_cm = plot_nm(X_test, y_test, model, br_title)
            fig_cm.savefig(f'{plot_dir}/cm_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_cm)
            
            # Plot ROC
            score_list, var_list, var_str = event_performance(all_df, model)
            fig_var_score = plot_var_score(var_list, score_list, var_str, f"Mass and Score (test, {br_title})")
            fig_var_score.savefig(f'{plot_dir}/pi0_mass_score_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_var_score)
            
            fig_roc = plot_roc(score_list, rf'ROC Curve - $\pi^{0}$ Classifier (test, {br_title})')
            fig_roc.savefig(f'{plot_dir}/roc_curv_{data_type}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_roc)
            
            print(f"\n✓ All results saved to {plot_dir}")
        
        else:
            print(f"Skipping {data_type} - not target category")
            