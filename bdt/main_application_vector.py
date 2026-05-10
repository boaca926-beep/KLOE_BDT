import pandas as pd
import numpy as np
import os
import sys
import joblib
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import seaborn as sns

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DATA_DIR, PLOT_APP_DIR, MODEL_DIR
from plots import plot_nm, plot_roc

def plot_event_confusion_matrix(event_results, data_type, plot_dir):
    cm = confusion_matrix(event_results['true_signal'], event_results['pred_signal'])
    cm_percent = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0],
                xticklabels=['Background', 'Signal'], yticklabels=['Background', 'Signal'])
    axes[0].set_title(f'Confusion Matrix (Counts) - {data_type}')
    sns.heatmap(cm_percent, annot=True, fmt='.1f', cmap='Blues', ax=axes[1],
                xticklabels=['Background', 'Signal'], yticklabels=['Background', 'Signal'])
    axes[1].set_title(f'Confusion Matrix (Percentages) - {data_type}')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/event_cm_{data_type}.png', dpi=300)
    plt.close(fig)
    tn, fp, fn, tp = cm.ravel()
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
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
    print("Application on test dataset (fast, exact original logic)...")
    input_data_dir = DATA_DIR
    input_model_dir = MODEL_DIR
    plot_dir = PLOT_APP_DIR
    import shutil
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir, exist_ok=True)

    phys_map = joblib.load(os.path.join(input_data_dir, 'phys_map.pkl'))
    category_type = 'TCOMB'

    for data_type, info in phys_map.items():
        if data_type != category_type:
            continue
        br_title = info['br_title']
        print(f"\nProcessing {data_type}...")

        all_df_test = joblib.load(os.path.join(input_data_dir, f'all_df_test_{data_type}.pkl'))
        X_test = joblib.load(os.path.join(input_data_dir, f'X_test_{data_type}.pkl'))
        y_test = joblib.load(os.path.join(input_data_dir, f'y_test_{data_type}.pkl'))
        model = joblib.load(os.path.join(input_model_dir, f'pi0_classifier_model_{data_type}.pkl'))

        print(f"Loaded: {len(all_df_test)} photons, {len(X_test)} pairs")

        y_pred_proba = model.predict_proba(X_test)[:, 1]
        print(f"Probabilities shape: {y_pred_proba.shape}")

        n_photons = len(all_df_test)
        n_pairs_per_photon = 3
        probs_photon = y_pred_proba.reshape(n_photons, n_pairs_per_photon)
        true_pair_labels = np.array(y_test).reshape(n_photons, n_pairs_per_photon)
        true_photon = (true_pair_labels.sum(axis=1) > 0).astype(int)   # 1 if photon has a true π⁰ pair

        max_proba_photon = probs_photon.max(axis=1)
        mean_proba_photon = probs_photon.mean(axis=1)
        event_ids = all_df_test['event'].values

        thresholds = np.arange(0.05, 1.0, 0.05)
        best_f1 = 0.0
        best_thr = 0.5
        best_strat = 'any'
        best_event_df = None

        for thr in thresholds:
            pred_binary = (probs_photon >= thr).astype(int)
            photon_has_any = (pred_binary.sum(axis=1) > 0).astype(int)

            photon_df = pd.DataFrame({
                'event_id': event_ids,
                'true_photon': true_photon,
                'photon_has_any': photon_has_any,
                'max_proba': max_proba_photon,
                'mean_proba': mean_proba_photon,
            })

            # Event aggregation: true_signal = OR over photons, pred_* from sums
            event_df = photon_df.groupby('event_id').agg({
                'true_photon': 'max',           # event true if any photon true
                'photon_has_any': 'sum',        # count of photons with any predicted pair
                'max_proba': 'max',
                'mean_proba': 'mean',
            }).reset_index()
            event_df.rename(columns={'true_photon': 'true_signal'}, inplace=True)

            event_df['pred_any'] = (event_df['photon_has_any'] > 0).astype(int)
            event_df['pred_min2'] = (event_df['photon_has_any'] >= 2).astype(int)
            event_df['pred_max'] = (event_df['max_proba'] >= thr).astype(int)
            event_df['pred_mean'] = (event_df['mean_proba'] >= thr).astype(int)

            best_f1_this = 0
            best_strat_this = 'any'
            for strat, col in [('any','pred_any'), ('min2','pred_min2'), ('max','pred_max'), ('mean','pred_mean')]:
                tp = ((event_df['true_signal'] == 1) & (event_df[col] == 1)).sum()
                fp = ((event_df['true_signal'] == 0) & (event_df[col] == 1)).sum()
                fn = ((event_df['true_signal'] == 1) & (event_df[col] == 0)).sum()
                precision = tp / (tp + fp) if (tp + fp) > 0 else 0
                recall = tp / (tp + fn) if (tp + fn) > 0 else 0
                f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
                if f1 > best_f1_this:
                    best_f1_this = f1
                    best_strat_this = strat

            if best_f1_this > best_f1:
                best_f1 = best_f1_this
                best_thr = thr
                best_strat = best_strat_this
                best_event_df = event_df.copy()
                best_event_df['pred_signal'] = best_event_df[f'pred_{best_strat_this}']

        print(f"\n{'='*60}")
        print("OPTIMAL CONFIGURATION")
        print(f"{'='*60}")
        print(f"Threshold  : {best_thr:.2f}")
        print(f"Strategy   : '{best_strat}'")
        print(f"F1 score   : {best_f1:.4f}")

        best_event_df.to_csv(f'{plot_dir}/event_results_{data_type}.csv', index=False)
        plot_event_confusion_matrix(best_event_df, data_type, plot_dir)

        fig_cm = plot_nm(X_test, y_test, model, br_title)
        fig_cm.savefig(f'{plot_dir}/cm_{data_type}.png')
        plt.close(fig_cm)

        fig_roc = plot_roc(y_test, y_pred_proba, rf'ROC Curve - π⁰ Classifier (test, {br_title})')
        fig_roc.savefig(f'{plot_dir}/roc_curv_{data_type}.png')
        plt.close(fig_roc)

        print(f"\n✓ All results saved to {plot_dir}")