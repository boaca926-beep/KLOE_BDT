# Validation metrics
import os
import sys
import joblib
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve, auc, accuracy_score
from prediction import find_best_pi0_candidate

'''
Your code is perfectly consistent! The 100% recall means the model always ranks the true π⁰ pair highest, but sometimes with low confidence (< threshold). 
The 90.7% recall is the rate at which the model's confidence exceeds threshold. Both metrics are valuable for understanding your model's performance!
'''

#================================================
# Evaluate validation set
#================================================
def eval_performance(model, X_val, y_val):
    y_pred_proba = model.predict_proba(X_val)[:, 1]
    y_pred = model.predict(X_val)

    auc = roc_auc_score(y_val, y_pred_proba)
    acc = accuracy_score(y_val, y_pred)

    print(f"\nValidation Performance:")
    print(f"AUC: {auc:.3f}")
    print(f"Accuracy: {acc:.4f}")

    return auc, acc

#================================================
# Evaluating performance
#================================================
def event_performance(all_df, model):
    print("\n* Evaluating performance on validation set:")

    result_events = all_df

    correct_predictions = 0
    total_pos_events = 0

    results = [] # Results list

    # Collect masses of identified pi0 candidates
    candidate_masses = []
    true_matches = []

    # Collect scores of identified pi0 candidates
    score_pos= []
    score_neg = []

    for _, evt in result_events.iterrows():
        # evt is a list object
        photons = [
            np.array([evt.Br_E1, evt.Br_px1, evt.Br_py1, evt.Br_pz1]),
            np.array([evt.Br_E2, evt.Br_px2, evt.Br_py2, evt.Br_pz2]),
            np.array([evt.Br_E3, evt.Br_px3, evt.Br_py3, evt.Br_pz3]),
        ]

        best_pair, score, mass = find_best_pi0_candidate(photons, model)
        
        # Create event_info dictionary for this event
        event_info ={
            'mass': mass,
            'score': score
        }

        # Save results
        results.append(event_info)

        # Couters
        if evt.is_signal:
            total_pos_events += 1
            #if best_pair == evt.true_pi0_pair:
            if np.array_equal(best_pair, evt.true_pi0_pair):

                correct_predictions += 1
                score_pos.append(score)
                status = "✓"
            else:
                score_neg.append(score)
                status = "✗"
        else:
            score_neg.append(score)
            status = "BG"

        truth = f"Postive pi0: {evt.true_pi0_pair}" if evt.is_signal else "Negative event"
        #print(f"Event {evt.event}: Best pair {best_pair}, score={score:.3f}, m={mass:.3f} | {status}")

        # mass collection
        if evt.is_signal and np.array_equal(best_pair, evt.true_pi0_pair):
        #if evt.is_signal and best_pair == evt.true_pi0_pair:
            candidate_masses.append(mass)
            true_matches.append(1)
        elif evt.is_signal:
            candidate_masses.append(mass)
            true_matches.append(0)

    if total_pos_events > 0:
        accuracy = correct_predictions / total_pos_events * 100 
        print(f"\nAccuracy on positive events: {accuracy:.1f}% ({correct_predictions}/{total_pos_events})")

    ## Check results list
    print(f"Collected {len(results)} events total")
    results_df = pd.DataFrame(results) # Covert to data frame
    print(results_df.head())

    ## Create score_pos and score_neg structure
    #print(rf'{len(score_pos)},     {len(score_neg)}')
    score_pos_df = pd.DataFrame({
        'score': score_pos,
        'type': 'signal_correct'
    })

    score_neg_df = pd.DataFrame({
        'score': score_neg,
        'type': 'background_or_wrong'
    })

    score_all_df = pd.concat([score_pos_df, score_neg_df], ignore_index=True) # Combining data from score_pos_df and score_neg_df
    #print(score_pos_df.head(5))
    #print(score_neg_df.head(5))
    #print(score_all_df.head(5))

    m_pos = [m for m, match in zip(candidate_masses, true_matches) if match]
    m_neg = [m for m, match in zip(candidate_masses, true_matches) if not match]
    #print(m_pos)
    print(m_neg)
    var_list = [m_pos, m_neg]
    score_list = [score_pos, score_neg]
    var_str = [rf'$\pi^{0}$ Candidate Mass Distribution', 'Events', rf'$M_{{\gamma\gamma}}$ ($\mathrm{{MeV}}/\mathrm{{c}}^{2}$)'] # [title, y_label, x_label]

    return score_list, var_list, var_str
    