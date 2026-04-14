# Plotting functions
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, auc, accuracy_score
from sklearn.metrics import confusion_matrix, classification_report

# =================================================================
# Plot (all, positive, negative) comparison
# =================================================================
def plot_compr_hist(df_set, drop_columns, rows, bins, plot_title):
                   
                   
    all_df = df_set[0].drop(drop_columns, axis=1)
    good_df = df_set[1].drop(drop_columns, axis=1)
    bad_df = df_set[2].drop(drop_columns, axis=1)

    ##  S/B ratio
    S = len(good_df)
    B = len(bad_df)
    S_purity = S / (S + B)
    print(f"Total events: {len(all_df)}, postive events: {len(good_df)}, negative events: {len(bad_df)}")
    print(f"S/sqrt(S+B): {S / np.sqrt(S + B):.2f}")
    print(f"S_purity: {S_purity:.2f}")

    ## Check col_len
    col_len = len(all_df.columns) # length of columns of df

    if (col_len < 0):
        # negative
        print(f"Negative col_len ({col_len})")
        return
    elif (col_len == 0):
        # zero col_len
        print(f"Zero col_len ({col_len})")
    else:
        # postive
        if (col_len % 2 == 0):
            # even case
            print(f"good events col_len ({col_len})")
        else:
            # odd or not integer
            print(all_df.columns)
            print(f"Odd column length or none integer column length ({col_len}). Not plot is created!")
            return

    # Create subplot grid
    plot_col = int(col_len / rows) # number of rows and columns to the plot
    fig, axes = plt.subplots(rows, plot_col, figsize=(16, 10)) # rows and columns to subplots
    fig.suptitle(plot_title, fontsize=16, y=1.02)

    # Flatten axes array for easy iteration
    axes = axes.flatten()
    columns = all_df.columns
    
    #for i, label in enumerate(columns_df[:col_len]):
    for i, label in enumerate(columns):
        #print(i, label)
        # desity=True normalized
        positive_good_df = good_df[label] #good_df[good_df[label] > 0.2][label]
        positive_bad_df = bad_df[label] #bad_df[bad_df[label] > 0.2][label]
        positive_all_df = all_df[label] #all_df[all_df[label] > 0.2][label]

        if label in ['Br_E1', 'Br_E2', 'Br_E3', 'Br_deltaE']:
            unit = fr'[$\mathrm{{MeV}}$]'
        elif label in ['Br_m_gg', 'Br_m3pi', 'Br_ppIM']:
            unit = fr'[$\mathrm{{MeV}}/\mathrm{{c}}^{2}$]'
        elif label in ['Br_px1', 'Br_py1', 'Br_pz1', 'Br_px2', 'Br_py2', 'Br_pz2', 'Br_px3', 'Br_py3', 'Br_pz3']:
            unit = fr'[$\mathrm{{MeV}}/\mathrm{{c}}$]'
        elif label in ['Br_angle_pi0gam12']:
            unit = fr'[$\circ$]'
        else:
            unit = ""
            #rint("AU")

        n1, bin_edges1, patches1 = axes[i].hist([positive_good_df, positive_bad_df], 
                     color=['green', 'blue'], 
                     bins=bins, 
                     label=[f'Positive', f'Negative'], 
                     density=False, 
                     edgecolor=['green', 'blue'],
                     linewidth=1, 
                     alpha=0.5,
                     histtype='stepfilled' # Filled histograms
                     )

        n2, bin_edges2, patches2 = axes[i].hist(positive_all_df, 
                     color=['red'], 
                     bins=bin_edges1, 
                     label='All', 
                     density=False, 
                     edgecolor='red',
                     linewidth=1, 
                     alpha=0.5,
                     histtype='step' # Filled histograms
                     )
        
        

        bin_width = bin_edges2[1] - bin_edges2[0]
        #print(f"bin_width: {bin_width:.2f}")
 
        axes[i].set_xlabel(label + ' ' + unit)
        #axes[i].set_ylabel(fr'Events / {bin_width:.1f} {unit}', fontsize=14)
        axes[i].set_ylabel(fr'Events', fontsize=14)
        axes[i].grid(True, alpha=0.3)
        axes[i].legend(loc='best', fontsize=14) 
    
    #plt.title(plot_title)
    plt.tight_layout()
    #plt.savefig('./plots/' + plot_nm + '_compr.png', dpi=300, bbox_inches='tight')
    plt.show(block=False)
    #plt.show()
    #plt.close()

    return fig

# =================================================================
# Plot single variable
# =================================================================
def plot_var(array, var_nm, phys_ch):
    print(f"Plotting ... {var_nm}")
    fig, ax = plt.subplots(figsize=(16, 10))
    plt.hist(array, 
             color='green', 
             bins=400, 
             density=False, 
             edgecolor='black', 
             alpha=0.7, 
             label=r'True postive',
             histtype='stepfilled'
             )
    plt.xlabel(r'$M_{\gamma\gamma}$ $[MeV/c^{2}]$', fontsize=14)                                  
    plt.ylabel('Events', fontsize=14)
    plt.title(fr'Mass Distribution of $M_{{\gamma\gamma}}$ (n={len(array)}) {phys_ch}', fontsize=16)
    # combine into one legend
    #plt.legend(loc='best', fontsize=8, frameon=True, fancybox=True, shadow=True,
    #           title=f'π⁰ Mass Distribution (n={len(array)})\nTrue π⁰ events'
    #)
    #plt.legend(loc='best', fontsize=8, title=f'π⁰ Mass Distribution (n={len(array)})') 
    plt.legend(loc='best', fontsize=14, frameon=True, fancybox=True, shadow=True)
    plt.grid(True, alpha=0.3)
    #plt.savefig('./plots/signal_pi0.png')
    plt.show(block=False)
    #plt.show()
    #plt.close()

    return fig

# =================================================================
# Plot feature-feature
# =================================================================
def plot_feature_pairs(df, drop_columns, plot_title, hue_tmp):
    print('Plotting feature pairs')

    feature_columns = [col for col in df.columns if col not in drop_columns]
    
    g = sns.pairplot(df[feature_columns], # Data
                     hue = hue_tmp, # Color grouping, points by the values in the 'is_pi0' column
                     palette={1: 'blue', 0: 'red'}, # 3. colors     
                     diag_kind='hist', # Diagonal plot type
                     plot_kws={'alpha': 0.5, 's': 10}, # Scatter plot options
                     diag_kws={'alpha': 0.7, 'edgecolor': 'black'} # Histogram options  
    )
    g.figure.suptitle(plot_title, y=1.02, fontsize=14)
    plt.tight_layout()
    #plt.savefig('./plots/' + plot_nm, dpi=300, bbox_inches='tight')
    plt.show(block=False)
    #plt.show()
    #plt.close()

    return g

# =================================================================
# Plot feature-target    
# ================================================================= 
def plot_feature_target(target_corr, plot_title):
    print('here plotting ...')
    fig, ax = plt.subplots(figsize=(10, 6))
    #plt.figure(figsize=(10, 6))
    target_corr_pos = [np.abs(e) for e in target_corr.values] # abs values
    colors = ['red']
    #colors = ['green' if c > 0  else 'red' for c in target_corr.values]
    #plt.bar(range(len(target_corr_pos)), target_corr_pos.values, color=colors, alpha=0.7)
    plt.bar(range(len(target_corr_pos)), target_corr_pos, color=colors, alpha=0.7)
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    plt.xticks(range(len(target_corr_pos)), target_corr.index, rotation=45, ha='right', fontsize=14)
    plt.ylabel(rf'Absolute correlation with true $\pi^{0}$', fontsize=14)
    plt.title(plot_title, fontsize=14)
    #plt.title(rf'Feature Importance: Correlation with true $\pi^{0}$')
    plt.grid(True, alpha=0.3, axis='y')

    # Add value labels
    for i, (feat, corr) in enumerate(target_corr.items()):
        #print (i, feat, corr)
        corr = np.abs(corr)
        plt.text(i, corr + (0.02 if corr > 0 else -0.05),
                    f'{corr:.2f}', ha='center', va='bottom' if corr > 0 else 'top')
        
    plt.tight_layout()
    #plt.savefig('./plots/' + plot_nm + '.png', dpi=300, bbox_inches='tight')
    plt.show(block=False)
    #plt.close()

    return fig


# =================================================================
# Plot variable vs. score 
# ================================================================= 
def plot_var_score(var_list, score_list, var_str, plot_title):
    print("Plotting variable vs. score ...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    fig.suptitle(plot_title, fontsize=16, y=1.02)

    titles = [var_str[0], 'Score']
    y_labels = [var_str[1], 'Events']
    x_labels = [var_str[2], 'Score']

    for i in range(2):
        if (i == 0): # mass distributions
            n, bin_edges, patches = axes[i].hist(var_list,
                                                bins=200, 
                                                alpha=0.5, 
                                                label=['Correctly identified', 'Wrongly identified'],
                                                color=['green', 'black'],
                                                density=False,
                                                linewidth=1,
                                                histtype='stepfilled'
                                                )
            
            axes[i].set_title(fr'{titles[i]}', fontsize=18)
            #axes[i].set_xlim(50, 200) # Set x-axis range in [MeV/c^2]
            axes[i].set_xlabel(fr'{x_labels[i]}', fontsize=14)
            axes[i].set_ylabel(fr'{y_labels[i]}', fontsize=14)
            axes[i].grid(True, alpha=0.3)
            axes[i].legend(loc='best', fontsize=14)
            #axes[i].set_yscale('log')
            axes[i].axvline(x=135, color='black', linestyle='--', label='True pi0 mass')
        else:
            n, bin_edges, patches = axes[i].hist(score_list,
                                                bins=100, 
                                                alpha=0.5, 
                                                label=['Signal', 'Background'],
                                                color=['blue', 'red'],
                                                density=False,
                                                linewidth=1,
                                                histtype='step'
                                                )
            axes[i].set_title(fr'{titles[i]}', fontsize=18)
            #axes[i].set_xlim(0, 0.2) # Set x-axis range from 0 to 0.2
            axes[i].set_xlabel(fr'{x_labels[i]}', fontsize=14)
            axes[i].set_ylabel(fr'{y_labels[i]}', fontsize=14)
            axes[i].grid(True, alpha=0.3)
            axes[i].legend(loc='best', fontsize=14)
            axes[i].axvline(x=0.5, color='black', linestyle='--', label='True pi0 mass')
            axes[i].set_yscale('log')

    plt.tight_layout()
    #plt.savefig(rf'./plots/{plot_nm}.png')
    plt.show(block=False)
    #plt.close()

    return fig

# =================================================================
# Plot ROC curve (Performance)
# =================================================================
def plot_roc(score_list, plot_title):
    print("Plotting ROC curv ...")
    """
    - score_pos: scores from positive class (correctly identified pi0)
    - score_neg: scores from negative class (wrongly identified pi0)
    """
    score_pos = score_list[0]
    score_neg = score_list[1]
    y_true = [1] * len(score_pos) + [0] * len(score_neg)
    y_score = score_pos + score_neg
    #print(len(y_true), y_true)
    #print(y_score)

    # Calculate ROC curve
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(10, 8))

    plt.plot(fpr,
             tpr,
             color='black',
             lw=2,
             label=f'ROC AUC: {roc_auc:.4f}'
             )
    plt.plot([0, 1], [0, 1], 'r--', label='Random')
    plt.xlabel('False Positive Rate', fontsize=14)
    plt.ylabel('True Positive Rate', fontsize=14)
    plt.title(plot_title)
    plt.legend(loc='best', fontsize=14)
    plt.grid(True, alpha=0.3)

    # Add statistics
    textstr = f'Positive: {len(score_pos)} events\nBackground: {len(score_neg)} events'
    plt.text(0.6, 0.2, textstr, fontsize=14,
             bbox=dict(boxstyle="round,pad=0.5",
             facecolor='yellow',
             alpha=0.3))
    plt.tight_layout()
    #plt.savefig(f'./plots/{plot_nm}.png', dpi=300)
    plt.show(block=False)
    #plt.close()

    return fig

# =================================================================
# Plot learning curves (Check for overfitting)
# =================================================================
def plot_learning_curves(model, plot_title):
    """
    Plot training vs validation performance over boosting rounds
    This shows if model is overfitting
    """
    print("Plotting learning curves...")

    early_stop = model.get_params()['early_stopping_rounds']
    print(early_stop)

    results = model.evals_result()
    #print(results)

    train_auc = results['validation_0']['auc']
    val_auc = results['validation_1']['auc']

    # Error rate (convert to accuracy)
    train_error = results['validation_0']['error']
    val_error = results['validation_1']['error']
    train_acc = [1 - err for err in train_error]
    val_acc = [1 - err for err in val_error]
    #print(train_auc)

    fig, axes = plt.subplots(1, 1, figsize=(10, 8))

    # AUC over rounds
    axes.plot(train_auc, 'b-', label='Training', linewidth=2)
    axes.plot(val_auc, 'r-', label='Validation', linewidth=2)
    #axes.set_ylim(0.8, 1) # Set y-axis range from 0 to 1
    axes.set_xlabel('Boosting Round', fontsize=14)
    axes.set_ylabel('AUC', fontsize=14)
    #axes.set_title(rf'Learning Curves - AUC', fontsize=14)
    axes.axvline(x=early_stop, color='black', linestyle='--', linewidth=2, 
               label=f'Early stop (iteration {early_stop})')
    axes.legend()
    axes.grid(True, alpha=0.3)


    plt.suptitle(plot_title, fontsize=14)
    plt.tight_layout()
    #plt.savefig(f'./plots/{plot_nm}.png', dpi=300)
    plt.show(block=False)
    #plt.close()

    # Print diagnostics
    final_gap = train_auc[-1] - val_auc[-1]
    print(f"\n📊 Validation Diagnostics:")
    print(f"  Final Training AUC: {train_auc[-1]:.4f}; Accuracy: {train_acc[-1]:.4}")
    print(f"  Final Validation AUC: {val_auc[-1]:.4f}; Accuracy: {val_acc[-1]:.4}")
    print(f"  Gap: {final_gap:.4f}")
    
    if final_gap > 0.05:
        print("  ⚠️  WARNING: Possible overfitting!")
    elif final_gap > 0.02:
        print("  ⚠️  Caution: Moderate gap")
    else:
        print("  ✅ Good generalization!")
    
    return fig

# =================================================================
# Plot confusion matrix
# =================================================================
def plot_nm(X_test, y_test, model, phys_ch):
    print("Plotting confusion matrix ...")

    y_pred = model.predict(X_test) # Prediction
    cm = confusion_matrix(y_test, y_pred) # Confusion matrix
    #print(cm)
    #print(X_test.columns)

    # Visualize it
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
    plt.colorbar(im, ax=ax) # Add color bar
    
    ax.set_title(rf'Confusion Matrix (test, {phys_ch})', fontsize=16)
    ax.set_xlabel('Predicted', fontsize=14)
    ax.set_ylabel('True', fontsize=14)

    # For the counts plot (ax[0])
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(['Background', 'Signal'])
    ax.set_yticklabels(['Background', 'Signal'])

    # Add text annotations
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, str(cm[i, j]), ha='center', va='center', color='red')
    plt.show(block=False)

    # Get detailed metrics
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))

    return fig