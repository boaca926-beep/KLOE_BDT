def get_cm(nb_tot, tp, fp, tn, fn):
    '''
    Store confusion matrix entries
    '''

    nb_sum = tn + fp + fn + tp
    accuracy = (tn + tp) / nb_sum
    precision = tp / (fp + tp)
    recall = tp / (fn + tp)
    F1 = 2 * (precision * recall) / (precision + recall)

    

    cm = {  
        'nb_tot': nb_tot,
        'tp': tp,
        'fp': fp,
        'tn': tn,
        'fn': fn,
        'cm_entr_sum': nb_sum,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'F1': F1
    }
    
    if nb_sum != nb_tot:
        print("Total number of photons and the totla number of entries are not equal!")
        return cm
    

    return cm

## Event-wise
#nb_evnt_val = 220463 # total number of photons
#tp = 141385 # false-negative 
#fp = 31035 # false-positive
#tn = 11450 # false-negative
#fn = 16349 # true-positive

cm_val = get_cm(220463, 141385, 31035, 11450, 16349) # tp, fp, tn, fn, event-wise
#cm_val = get_cm(661389, 148905, 42145, 449087, 21252) # tp, fp, tn, fn, pair-wise, theshold=0.5

print(type(cm_val))
#fp = cm_evnt_

tn_val=cm_val['tn']
fp_val=cm_val['fp']
fn_val=cm_val['fn']
tp_val=cm_val['tp']

cm_entr_sum_val = cm_val['cm_entr_sum'] #tn + fp + fn + tp
accuracy = cm_val['accuracy'] #(tn + tp) / cm
precision = cm_val['precision'] #tp / (fp + tp)
recall = cm_val['recall'] #tp / (fn + tp)
F1 = 2 * (precision * recall) / (precision + recall)

print(f"Total number of validation events: {cm_entr_sum_val}\n")
print(f"    1: {tp_val + fn_val} ({(tp_val + fn_val)/cm_entr_sum_val*100:0.1f}%)")
print(f"    0: {tn_val + fp_val} ({(tn_val + fp_val)/cm_entr_sum_val*100:0.1f}%)\n")

#print(f"    1+0: {nb_evnt_val_pair_bad + nb_evnt_val_pair_good}\n")
print(f"    tn: {tn_val} fp: {fp_val}")
print(f"    fn: {fn_val} tp: {tp_val}\n")
print(f"    cm entry: {cm_entr_sum_val}")
print(f"    accuracy = {accuracy* 100:0.1f}%")
print(f"    precision = {precision* 100:0.1f}%")
print(f"    recall = {recall* 100:0.1f}%")
print(f"    F1 = {F1* 100:0.1f}%") # Harmonic Mean = 2 / (1/a + 1/b) = 2ab / (a + b)


r'''
# Pair-wise
nb_evnt_val_pair = 661389 # total number of test pairs
nb_evnt_val_pair_good = 170157 # true signal events (kloe selection)
nb_evnt_val_pair_bad = 491232 # background events (kloe selection)

## Confusion matrix
tn_pair = 449087 # false-negative 
fp_pair = 42145 # false-positive
fn_pair = 21252 # false-negative
tp_pair = 148905 # true-positive
cm_pair_entry = tn_pair + fp_pair + fn_pair + tp_pair


print(f"Total number of validation events: {nb_evnt_val}, pairs: {nb_evnt_val_pair}\n")
print(f"Number of pairs per events: {nb_evnt_val_pair/nb_evnt_val}")
print(f"    1: {nb_evnt_val_pair_good} ({nb_evnt_val_pair_good/nb_evnt_val_pair*100:0.1f}%) 0: {nb_evnt_val_pair_bad} ({nb_evnt_val_pair_bad/nb_evnt_val_pair*100:0.1f}%)")
print(f"    1+0: {nb_evnt_val_pair_bad + nb_evnt_val_pair_good}\n")
print(f"    tn: {tn_pair} fp: {fp_pair}")
print(f"    fn: {fn_pair} tp: {tp_pair}")
print(f"cm entry: {cm_pair_entry}")
'''