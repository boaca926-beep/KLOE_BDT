# Plotting functions
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, auc, accuracy_score
from sklearn.metrics import confusion_matrix, classification_report

#uv run main_validation.py 2>&1 | tee validation_log.txt

# =================================================================
# Plot (all, positive, negative) comparison
# =================================================================
def plot_compr_hist(df_set, drop_columns, rows=3, bins=50, plot_title="", 
                    subplot_size=4, units=None):
    """
    Plot histograms comparing Signal vs Background for all features.
    
    Args:
        df_set: list of [all_df, signal_df, background_df]
        drop_columns: columns to exclude from plotting
        rows: number of columns per row (default: 3)
        bins: number of histogram bins
        plot_title: title for the figure
        subplot_size: size in inches for each subplot (width = height)
        units: optional dict for custom units
    """
    all_df = df_set[0].drop(drop_columns, axis=1)
    good_df = df_set[1].drop(drop_columns, axis=1)
    bad_df = df_set[2].drop(drop_columns, axis=1)

    ## S/B ratio
    S = len(good_df)
    B = len(bad_df)
    S_purity = S / (S + B)
    print(f"Total events: {len(all_df)}, signal: {len(good_df)}, background: {len(bad_df)}")
    print(f"S/sqrt(S+B): {S / np.sqrt(S + B):.2f}")
    print(f"S_purity: {S_purity:.2f}")

    col_len = len(all_df.columns)

    # ========== Define display names and units ==========
    display_names = {
        'm_gg': r'$m_{\gamma\gamma}$',
        'opening_angle': r'$\angle_{\gamma\gamma}$',
        'cos_theta': r'$\cos\theta_{\gamma\gamma}$',
        'E_asym': r'$A_{E}$',
        'e_min_x_angle': r'$\min(E_{\gamma_{1}},E_{\gamma_{2}})\times\theta_{\gamma\gamma}$',
        'E1': r'$E_{\gamma_{1}}$',
        'E2': r'$E_{\gamma_{2}}$',
        'E3': r'$E_{\gamma_{3}}$',
        'E_diff': r'$\left|E_{\gamma_{1}}-E_{\gamma_{2}}\right|$',
        'asym_x_angle': r'$A_{E}\times\theta_{\gamma\gamma}$',
        'Br_E1': r'$E_{\gamma_{1}}$',
        'Br_E2': r'$E_{\gamma_{2}}$',
        'Br_E3': r'$E_{\gamma_{3}}$',
        'Br_px1': r'$p_{x_{1}}$',
        'Br_py1': r'$p_{y_{1}}$',
        'Br_pz1': r'$p_{z_{1}}$',
        'Br_px2': r'$p_{x_{2}}$',
        'Br_py2': r'$p_{y_{2}}$',
        'Br_pz2': r'$p_{z_{2}}$',
        'Br_px3': r'$p_{x_{3}}$',
        'Br_py3': r'$p_{y_{3}}$',
        'Br_pz3': r'$p_{z_{3}}$',
        'Br_Eprompt_max': r'$E^{max}_{\gamma}$',
        'Br_m3pi': r'$M_{3\pi}$',
        'Br_ppIM': r'$M_{\pi^+\pi^-}$',
        'Br_deltaE': r'$\Delta E$',
        'Br_angle_pi0gam12': r'$\theta_{\pi^0\gamma}$',
        'Br_betapi0': r'$\beta_{\pi^0}$',
        'Br_lagvalue_min_7C': r'$\chi^{2}_{7C}$',
    }

    unit_map = {
        'm_gg': r'[MeV/$c^2$]',
        'opening_angle': r'[rad]',
        'cos_theta': r'',
        'E_asym': r'',
        'e_min_x_angle': r'[MeV]',
        'asym_x_angle': r'[MeV]',
        'E1': r'[MeV]',
        'E2': r'[MeV]',
        'E3': r'[MeV]',
        'E_diff': r'[MeV]',
        'Br_E1': r'[MeV]',
        'Br_E2': r'[MeV]',
        'Br_E3': r'[MeV]',
        'Br_px1': r'[MeV/$c$]',
        'Br_py1': r'[MeV/$c$]',
        'Br_pz1': r'[MeV/$c$]',
        'Br_px2': r'[MeV/$c$]',
        'Br_py2': r'[MeV/$c$]',
        'Br_pz2': r'[MeV/$c$]',
        'Br_px3': r'[MeV/$c$]',
        'Br_py3': r'[MeV/$c$]',
        'Br_pz3': r'[MeV/$c$]',
        'Br_Eprompt_max': r'[MeV]',
        'Br_m3pi': r'[MeV/$c^2$]',
        'Br_ppIM': r'[MeV/$c^2$]',
        'Br_deltaE': r'[MeV]',
        'Br_angle_pi0gam12': r'[$^\circ$]',
        'Br_betapi0': r'',
        'Br_lagvalue_min_7C': r'',
    }
    
    if units:
        unit_map.update(units)

    # ========== Calculate grid dimensions ==========
    plot_col = rows
    plot_row = (col_len + plot_col - 1) // plot_col
    
    # ========== FIX: Dynamic figure size for square subplots ==========
    fig_width = subplot_size * plot_col
    fig_height = subplot_size * plot_row
    fig, axes = plt.subplots(plot_row, plot_col, figsize=(fig_width, fig_height))
    # ================================================================
    
    fig.suptitle(plot_title, fontsize=16, y=1.02)

    axes = axes.flatten()
    columns = all_df.columns
    
    for i in range(col_len):
        label = columns[i]
        
        signal_vals = good_df[label].dropna()
        bkg_vals = bad_df[label].dropna()
        all_vals = all_df[label].dropna()

        unit = unit_map.get(label, "")
        display_name = display_names.get(label, label.replace('Br_', ''))

        axes[i].hist([signal_vals, bkg_vals], 
                     color=['green', 'blue'], 
                     bins=bins, 
                     label=['Signal', 'Background'], 
                     density=False, 
                     edgecolor=['green', 'blue'],
                     linewidth=1, 
                     alpha=0.4,
                     histtype='stepfilled'
        )
        
        #axes[i].hist(all_vals, 
        #             color='red', 
        #             bins=bins, 
        #             label='All', 
        #             density=False, 
        #             edgecolor='red',
        #             linewidth=1.5, 
        #             alpha=0.0,
        #             histtype='step'
        #)

        # Add statistics
        n_signal = len(signal_vals)
        n_bkg = len(bkg_vals)
        #axes[i].text(0.02, 0.95, f'S={n_signal}\nB={n_bkg}', 
        #            transform=axes[i].transAxes, fontsize=8,
        #            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        if unit:
            xlabel = display_name + ' ' + unit
        else:
            xlabel = display_name
        axes[i].set_xlabel(xlabel, fontsize=12)
        axes[i].set_ylabel('Events', fontsize=12)
        
        if i == 0:
            axes[i].legend(loc='upper left', fontsize=10)
        
        axes[i].grid(True, alpha=0.2)

    # Hide unused subplots
    for i in range(col_len, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
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
    #plt.show(block=False)
    #plt.show()
    #plt.close()

    return fig

# =================================================================
# Plot feature-feature
# =================================================================
def plot_feature_pairs(df, drop_columns, plot_title, hue_tmp, 
                       sample_frac=1.0, random_state=42, log_columns=None):
    """
    Create a pairplot with improved visibility and custom labels.
    
    Args:
        df: DataFrame
        drop_columns: columns to exclude from plotting
        plot_title: figure title
        hue_tmp: column name for color grouping
        sample_frac: fraction of data to use (for large datasets)
        random_state: random seed for sampling
        log_columns: list of column names to apply log scale (off-diagonal)
    """
    print('Plotting feature pairs')

    # Sort so signal (1) is plotted on top of background (0)
    df = df.sort_values(by=hue_tmp, ascending=True)

    if sample_frac < 1.0:
        df = df.sample(frac=sample_frac, random_state=random_state)
        print(f"Sampled {len(df)} rows for plotting")

    # Exclude both drop_columns AND hue_tmp from features
    feature_columns = [col for col in df.columns if col not in drop_columns and col != hue_tmp]
    
    if not feature_columns:
        print("Warning: No feature columns to plot!")
        return None
    
    print("df columns: ", df.columns.tolist())
    print(f"Feature columns: {feature_columns}")
    print(f"Hue column: {hue_tmp}")

    # ========== Define display names and units ==========
    display_names = {
        'pi0_label': r'True $\pi^{0}$',          # for custom label column
        'is_pi0': r'True $\pi^{0}$',             # fallback
        'Br_betapi0': r'$\beta_{\pi^0}$',
        'Br_ppIM': r'$M_{\pi^+\pi^-}$',
        'Br_angle_pi0gam12': r'$\theta_{\pi^0\gamma}$',
        'Br_deltaE': r'$\Delta E$',
        'Br_Eprompt_max': r'$E^{max}_{\gamma}$',
        'Br_m3pi': r'$M_{3\pi}$',
        'Br_lagvalue_min_7C': r'$\chi^{2}_{7C}$',
        'm_gg': r'$m_{\gamma\gamma}$',
        'opening_angle': r'$\angle_{\gamma\gamma}$',
        'cos_theta': r'$\cos\theta_{\gamma\gamma}$',
        'E_asym': r'$A_{E}$',
        'e_min_x_angle': r'$\min(E_{\gamma_{1}},E_{\gamma_{2}})\times\theta_{\gamma\gamma}$',
        'E1': r'$E_{\gamma_{1}}$',
        'E2': r'$E_{\gamma_{2}}$',
        'E3': r'$E_{\gamma_{3}}$',
        'E_diff': r'$\left|E_{\gamma_{1}}-E_{\gamma_{2}}\right|$'
    }

    unit_map = {
        'Br_betapi0': r'',
        'Br_ppIM': r'[MeV/$c^2$]',
        'Br_angle_pi0gam12': r'[$^\circ$]',
        'Br_deltaE': r'[MeV]',
        'Br_Eprompt_max': r'[MeV]',
        'Br_m3pi': r'[MeV/$c^2$]',
        'Br_lagvalue_min_7C': r'',
        'm_gg': r'[MeV/$c^2$]',
        'opening_angle': r'[rad]',
        'cos_theta': r'',
        'E_asym': r'',
        'e_min_x_angle': r'[MeV]',
        'asym_x_angle': r'[MeV]',
        'E1': r'[MeV]',
        'E2': r'[MeV]',
        'E3': r'[MeV]',
        'E_diff': r'[MeV]'
    }

    def get_label(col):
        """Get formatted label with display name and unit"""
        display_name = display_names.get(col, col)
        unit = unit_map.get(col, '')
        return f'{display_name} {unit}'.strip() if unit else display_name

    # Determine the palette based on the actual values in hue column
    unique_vals = df[hue_tmp].unique()
    print(f"Unique values in '{hue_tmp}': {unique_vals}")
    
    # Create appropriate palette
    if set(unique_vals) == {0, 1} or set(unique_vals) == {0.0, 1.0}:
        palette = {1: 'blue', 0: 'red'}
    elif set(unique_vals) == {True, False}:
        palette = {True: 'blue', False: 'red'}
    else:
        # For other values, create a generic palette
        palette = {val: 'blue' if i == 0 else 'red' for i, val in enumerate(sorted(unique_vals))}

    # Create the pairplot – pass the full DataFrame and specify vars
    g = sns.pairplot(df,
                     vars=feature_columns,
                     hue=hue_tmp,
                     palette=palette,
                     diag_kind='hist',
                     plot_kws={'alpha': 0.3, 's': 5, 'edgecolor': 'none'},  # reduced overlap
                     diag_kws={'alpha': 0.7, 'edgecolor': 'black'}
    )
    
    # Apply custom labels with units to the outer axes
    for i, col in enumerate(feature_columns):
        label = get_label(col)
        g.axes[i, 0].set_ylabel(label)   # left column
        g.axes[-1, i].set_xlabel(label)  # bottom row

    # Apply log scales to off-diagonal plots if requested
    if log_columns is not None:
        for i, col_x in enumerate(feature_columns):
            for j, col_y in enumerate(feature_columns):
                if i != j:
                    ax = g.axes[i, j]
                    if col_x in log_columns:
                        ax.set_xscale('log')
                    if col_y in log_columns:
                        ax.set_yscale('log')

    # --- Set legend title from display_names ---
    legend_title = display_names.get(hue_tmp, hue_tmp)
    if hasattr(g, 'legend'):
        g.legend.set_title(legend_title)
    elif hasattr(g, '_legend'):
        g._legend.set_title(legend_title)

    g.figure.suptitle(plot_title, y=1.02, fontsize=14)
    plt.tight_layout()
    
    return g

def plot_feature_pairs_old(df, drop_columns, plot_title, hue_tmp):
    print('Plotting feature pairs')

    feature_columns = [col for col in df.columns if col not in drop_columns]

    print("df columns: ", df.columns.tolist())
    print(f"Feature columns: {feature_columns}")
    print(f"Hue column: {hue_tmp}")
    
    df = df.sort_values(by=hue_tmp, ascending=True)

    g = sns.pairplot(df[feature_columns],                          # Data
                     hue = hue_tmp,                                # Color grouping, points by the values in the 'is_pi0' column
                     palette={1: 'blue', 0: 'red'},                # colors     
                     diag_kind='hist',                             # Diagonal plot type
                     plot_kws={'alpha': 0.5, 's': 10},             # Scatter plot options
                     diag_kws={'alpha': 0.5, 'edgecolor': 'black'} # Histogram options  
    )
    g.figure.suptitle(plot_title, y=1.02, fontsize=14)
    plt.tight_layout()
    #plt.savefig('./plots/' + plot_nm, dpi=300, bbox_inches='tight')
    #plt.show(block=False)
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
    #plt.show(block=False)
    #plt.close()

    return fig

def plot_feature_target_h(target_corr, plot_title):
    
    # ========== Define display names and units ==========
    display_names = {
        'pi0_label': r'True $\pi^{0}$',          # for custom label column
        'is_pi0': r'True $\pi^{0}$',             # fallback
        'Br_betapi0': r'$\beta_{\pi^0}$',
        'Br_ppIM': r'$M_{\pi^+\pi^-}$',
        'Br_angle_pi0gam12': r'$\theta_{\pi^0\gamma}$',
        'Br_deltaE': r'$\Delta E$',
        'Br_Eprompt_max': r'$E^{max}_{\gamma}$',
        'Br_m3pi': r'$M_{3\pi}$',
        'Br_lagvalue_min_7C': r'$\chi^{2}_{7C}$',
        'm_gg': r'$m_{\gamma\gamma}$',
        'opening_angle': r'$\angle_{\gamma\gamma}$',
        'cos_theta': r'$\cos\theta_{\gamma\gamma}$',
        'E_asym': r'$A_{E}$',
        'e_min_x_angle': r'$\min(E_{\gamma_{1}},E_{\gamma_{2}})\times\theta_{\gamma\gamma}$',
        'asym_x_angle': r'$A_{E}\times\theta_{\gamma\gamma}$',
        'E1': r'$E_{\gamma_{1}}$',
        'E2': r'$E_{\gamma_{2}}$',
        'E3': r'$E_{\gamma_{3}}$',
        'E_diff': r'$\left|E_{\gamma_{1}}-E_{\gamma_{2}}\right|$'
    }

    unit_map = {
        'Br_betapi0': r'',
        'Br_ppIM': r'[MeV/$c^2$]',
        'Br_angle_pi0gam12': r'[$^\circ$]',
        'Br_deltaE': r'[MeV]',
        'Br_Eprompt_max': r'[MeV]',
        'Br_m3pi': r'[MeV/$c^2$]',
        'Br_lagvalue_min_7C': r'',
        'm_gg': r'[MeV/$c^2$]',
        'opening_angle': r'[rad]',
        'cos_theta': r'',
        'E_asym': r'',
        'e_min_x_angle': r'[MeV]',
        'asym_x_angle': r'[MeV]',
        'E1': r'[MeV]',
        'E2': r'[MeV]',
        'E3': r'[MeV]',
        'E_diff': r'[MeV]'
    }

    def get_label(col):
        """Get formatted label with display name and unit"""
        display_name = display_names.get(col, col)
        unit = unit_map.get(col, '')
        return f'{display_name} {unit}'.strip() if unit else display_name
    
    print('Plotting feature correlations (horizontal)...')

    # Convert to Series if dict
    if isinstance(target_corr, dict):
        target_corr = pd.Series(target_corr)
    
    # Use absolute values for bar lengths
    corr_abs = np.abs(target_corr)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(corr_abs) * 0.4)))
    
    # Create horizontal bars
    bars = ax.barh(range(len(corr_abs)), corr_abs.values, 
                   color='red', alpha=0.7)
    
    # Set y-ticks with feature names
    ax.set_yticks(range(len(corr_abs)))
    #ax.set_yticklabels(target_corr.index, fontsize=12)
    ax.set_yticklabels([get_label(col) for col in target_corr.index],
                       rotation=45, ha='right', fontsize=12)
    
    # Axis labels and title
    ax.set_xlabel(r'Absolute correlation with true $\pi^{0}$', fontsize=14)
    ax.set_title(plot_title, fontsize=14)
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_xlim(0, 1)
    
    # Add value labels at the end of each bar
    for i, (idx, corr) in enumerate(corr_abs.items()):
        ax.text(corr + 0.01, i, f'{corr:.2f}', 
                va='center', ha='left', fontsize=10)
    
    # Invert y-axis to show highest correlation at top
    ax.invert_yaxis()
    
    plt.tight_layout()
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
    #plt.show(block=False)
    #plt.close()

    return fig

# =================================================================
# Plot ROC curve (Performance)
# =================================================================
def plot_roc(y_true, y_score, plot_title):
    """
    Plot ROC curve from true labels and predicted scores (probabilities).
    y_true: array-like of true binary labels (0/1)
    y_score: array-like of predicted probabilities for class 1 (signal)
    """
    print("Plotting ROC curve...")
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(fpr, tpr, lw=2, label=f'ROC AUC: {roc_auc:.4f}')
    ax.plot([0, 1], [0, 1], 'r--', label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=14)
    ax.set_ylabel('True Positive Rate', fontsize=14)
    ax.set_title(plot_title)
    ax.legend(loc='best', fontsize=14)
    ax.grid(True, alpha=0.3)

    # Optional: add counts
    n_pos = sum(y_true)
    n_neg = len(y_true) - n_pos
    textstr = f'Positive: {n_pos} events\nBackground: {n_neg} events'
    ax.text(0.6, 0.2, textstr, fontsize=14,
            bbox=dict(boxstyle="round,pad=0.5", facecolor='yellow', alpha=0.3))
    plt.tight_layout()
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
    #plt.show(block=False)
    #plt.close()

    # Print diagnostics
    final_gap = train_auc[-1] - val_auc[-1]
    print(f"\n Validation Diagnostics:")
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

def plot_learning_curves_improved(model, phys_ch, final_val_auc=None, final_gap=None, plot_title=""):
    """
    Enhanced learning curves plot with:
      - Zoomed y‑axis (0.94–1.0) to reveal tiny gaps.
      - Separate subplot for Train‑Validation gap.
      - Automatic annotation of early‑stop round (if any) and final performance.
    """
    # Extract data from the model
    results = model.evals_result()
    train_auc = np.array(results['validation_0']['auc'])
    val_auc   = np.array(results['validation_1']['auc'])
    rounds = np.arange(len(train_auc))

    # Error rate (convert to accuracy)
    train_error = results['validation_0']['error']
    val_error = results['validation_1']['error']
    train_acc = [1 - err for err in train_error]
    val_acc = [1 - err for err in val_error]

    # Retrieve early‑stop setting
    early_stop = model.get_params().get('early_stopping_rounds', None)
    # If early_stop is set, mark the exact round where it would have stopped
    # (usually the best iteration, but we'll just use the parameter value)
    # For the annotation, we'll use the round number given in the log.
    # If you have the best iteration from the model, use: model.best_iteration
    best_iter = getattr(model, 'best_iteration', None)
    if best_iter is None:
        best_iter = early_stop  # fallback

    # Compute gap
    gap = train_auc - val_auc

    #axes.set_ylim(0.8, 1) # Set y-axis range from 0 to 1
    y_start = 0.995
    x_start = 200

    # Create figure with two subplots (2 rows, 1 column)
    fig = plt.figure(figsize=(12, 9), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[2.5, 1], hspace=0.15)
    ax1 = fig.add_subplot(gs[0])  # main AUC
    ax2 = fig.add_subplot(gs[1])  # gap

    # ---------- MAIN AUC PLOT ----------
    ax1.tick_params(axis='both', labelsize=18)   # or 16, 20, etc.
    ax2.tick_params(axis='both', labelsize=18)

    ax1.set_ylim(y_start, .999)  # zoom to show all variations
    ax1.set_ylabel('AUC Score', fontsize=18, fontweight='bold')
    ax1.set_title(plot_title, fontsize=18, fontweight='bold', pad=15)
    #ax1.set_title(plot_title or f'Learning Curve ({phys_ch})', fontsize=18, fontweight='bold', pad=15)

    ax1.plot(rounds, train_auc, label='Training AUC', color='#1f77b4', linewidth=2.5, marker='o', markersize=4)
    ax1.plot(rounds, val_auc,   label='Validation AUC', color='#ff7f0e', linewidth=2.5, marker='s', markersize=4)

    # Annotate early‑stop (if given)
    if early_stop is not None and early_stop < len(rounds):
        ax1.axvline(x=early_stop, color='red', linestyle='--', alpha=0.6, linewidth=1.5)
        ax1.text(early_stop, y_start, 'Early stop', color='red', ha='center', fontsize=16, fontweight='bold',
                 bbox=dict(facecolor='white', edgecolor='red', boxstyle='round,pad=0.3'))

    # Annotate final convergence point (last round)
    x_shift = 30
    last_round = rounds[-1]
    final_val = val_auc[-1]
    ax1.scatter([last_round], [final_val], color='green', s=120, zorder=5, marker='*', edgecolors='black')
    ax1.text(last_round - x_shift, final_val - 0.0005, f'Final AUC = {final_val:.4f}', 
             color='darkgreen', ha='center', fontsize=16, fontweight='bold')

    ax1.grid(True, linestyle=':', alpha=0.7)
    ax1.legend(loc='lower right', fontsize=16, framealpha=0.95)
    ax1.set_xticklabels([])  # hide x‑ticks on top plot

    # ---------- GAP SUBPLOT ----------
    ax2.set_ylim(-0.0005, 0.002)  # tight around expected gap
    ax2.set_ylabel('Train - Val Gap', fontsize=18, fontweight='bold')
    ax2.set_xlabel('Boosting Round', fontsize=18, fontweight='bold')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)

    ax2.fill_between(rounds, 0, gap, color='skyblue', alpha=0.4, label='Gap area')
    ax2.plot(rounds, gap, color='navy', linewidth=2, marker='d', markersize=4, label='Gap size')

    # Annotate final gap (if provided, else compute)
    final_gap_val = final_gap if final_gap is not None else gap[-1]
    ax2.scatter([last_round], [gap[-1]], color='green', s=100, zorder=5, marker='*')
    ax2.text(last_round - x_shift, gap[-1] + 0.0008, f'Gap = {gap[-1]:.4f}', 
             color='darkgreen', ha='center', fontsize=16, fontweight='bold')

    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.legend(loc='upper left', fontsize=16)

    # Adjust x‑limits to give padding
    x_pad = max(2, int(0.03 * x_start))  # ~3% of 200 = 6

    ax1.set_xlim(-x_pad, len(rounds) + x_pad)
    ax2.set_xlim(-x_pad, len(rounds) + x_pad)
    
    # Print diagnostics
    final_gap = train_auc[-1] - val_auc[-1]
    print(f"\n Validation Diagnostics:")
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
def plot_cm(X_test, y_test, model, plot_title=""):
    print("Plotting confusion matrix ...")

    y_pred = model.predict(X_test) # Prediction
    cm = confusion_matrix(y_test, y_pred) # Confusion matrix
    #print(cm)
    #print(X_test.columns)

    # Visualize it
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
    plt.colorbar(im, ax=ax) # Add color bar
    
    #ax.set_title(rf'Confusion Matrix (test, {phys_ch})', fontsize=16)
    ax.set_xlabel('Predicted', fontsize=14)
    ax.set_ylabel('True', fontsize=14)
    ax.set_title(plot_title, fontsize=14, fontweight='bold', pad=15)

    # For the counts plot (ax[0])
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(['Background', 'Signal'])
    ax.set_yticklabels(['Background', 'Signal'])

    # Add text annotations
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, str(cm[i, j]), ha='center', va='center', color='red')
    #plt.show(block=False)

    # Get detailed metrics
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))

    return fig

def plot_cm_improved(X_test, y_test, model, phys_ch, plot_title=""):
    """
    Presentation-ready confusion matrix with:
      - Raw counts + row-normalized percentages in each cell.
      - Dynamic text color (white on dark cells, black on light cells).
      - No colorbar (saves space).
      - Larger, bold fonts.
    """
    print("Plotting confusion matrix (presentation layout)...")
    
    y_pred = model.predict(X_test)
    cm = confusion_matrix(y_test, y_pred)
    
    # Row-normalized percentages (efficiency per true class)
    cm_norm = cm.astype('float') / cm.sum(axis=1, keepdims=True)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Use the raw counts for the colormap intensity
    im = ax.imshow(cm, cmap='Blues', interpolation='nearest')
    
    # --- No colorbar for 2x2 (comment out if you really want it) ---
    # plt.colorbar(im, ax=ax)
    
    # Axis labels
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])


    ax.set_xticklabels(['Background', 'Signal'], fontsize=14)
    ax.set_yticklabels(['Background', 'Signal'], fontsize=14, rotation=90)
    ax.set_xlabel('Predicted', fontsize=16, fontweight='bold')
    ax.set_ylabel('True', fontsize=16, fontweight='bold')
    ax.set_title(plot_title or f'Confusion Matrix ({phys_ch})', fontsize=18, fontweight='bold', pad=20)
    
    # --- Add text with counts + percentages ---
    for i in range(2):
        for j in range(2):
            count = cm[i, j]
            pct = cm_norm[i, j] * 100
            
            # Choose text color based on cell intensity (dynamic contrast)
            # If the count is > 50% of the row total, use white; else black.
            row_total = cm[i, :].sum()
            threshold = row_total * 0.5
            text_color = 'white' if count > threshold else 'black'
            
            # Format the label: count on first line, percentage on second
            label = f"{count:,}\n({pct:.1f}%)"
            
            ax.text(j, i, label, ha='center', va='center', 
                    color=text_color, fontsize=14, fontweight='bold')
    
    # Optional: Add a grid line between cells for clarity
    ax.grid(False)
    # Draw lines to separate the 4 quadrants
    ax.axhline(y=0.5, color='black', linewidth=1.5)
    ax.axhline(y=1.5, color='black', linewidth=1.5)
    ax.axvline(x=0.5, color='black', linewidth=1.5)
    ax.axvline(x=1.5, color='black', linewidth=1.5)
    
    plt.tight_layout()
    return fig

def plot_event_cm(event_results, data_type):
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
    
    plt.tight_layout()
    return fig

def plot_event_cm_improved(event_results, data_type):
    print("Plotting event confusion matrix (presentation layout)...")
    y_true = event_results['true_signal']
    y_pred = event_results['pred_signal']
    cm = confusion_matrix(y_true, y_pred)
    cm_norm = cm.astype('float') / cm.sum(axis=1, keepdims=True)

    # Calculate metrics
    tn, fp, fn, tp = cm.ravel()
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor('#e8e8e8')   # light grey background

    im = ax.imshow(cm, cmap='Blues', interpolation='nearest',
                   vmin=0, vmax=cm.max() * 1.0)   # no extra padding

    # --- CORRECTED: Set ticks at cell centers ---
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(['Background', 'Signal'], fontsize=16)
    ax.set_yticklabels(['Background', 'Signal'], fontsize=16, rotation=90)
    
    # --- Ensure labels are centered on cells ---
    # This is the key: align labels to the center of the cells
    ax.set_xticks([0, 1], minor=False)
    ax.set_yticks([0, 1], minor=False)
    
    ax.set_xlabel('Predicted', fontsize=20, fontweight='bold')
    ax.set_ylabel('True', fontsize=20, fontweight='bold')

    max_count = cm.max()
    for i in range(2):
        for j in range(2):
            count = cm[i, j]
            pct = cm_norm[i, j] * 100
            text_color = 'white' if count > 0.2 * max_count else 'black'
            label = f"{count:,}\n({pct:.1f}%)"
            ax.text(j, i, label, ha='center', va='center',
                    color=text_color, fontsize=16, fontweight='bold')

    # Grid lines
    ax.axhline(y=0.5, color='black', linewidth=1.5)
    ax.axhline(y=1.5, color='black', linewidth=1.5)
    ax.axvline(x=0.5, color='black', linewidth=1.5)
    ax.axvline(x=1.5, color='black', linewidth=1.5)
    ax.grid(False)

    plt.tight_layout()
    return fig, accuracy, precision, recall, f1

# =================================================================
# Plot mass signal breakdown
# =================================================================
def plot_mass_signal_breakdown(X_test, y_test, y_pred, phys_ch="", plot_title=""):
    """
    Plot m_gg distribution split into:
        1. True Positives (Correct Signal)
        2. False Negatives (Signal Misclassified)
        3. True Negatives (Correct Background)
        4. False Positives (Background Misclassified)
        5. All physical background (True Negatives + False Positives)
    
    Args:
        X_test: DataFrame of test features (must contain 'm_gg')
        y_test: True labels (1 = signal, 0 = background)
        y_pred: Predicted labels
        phys_ch: Physics channel name (for title)
        plot_title: Custom title (overrides default)
    """
    print("Plotting mass distribution: good signal vs lost signal vs background...")
    
    # --- Split the data ---
    mask_tp = (y_test == 1) & (y_pred == 1)   # Good signal (True Positives)
    mask_fn = (y_test == 1) & (y_pred == 0)   # Lost signal (False Negatives)
    mask_bkg = (y_test == 0)                  # All physical background (True Negatives + False Positives)
    mask_sig = (y_test == 1)                  # All background (TP + FN)

    mass_tp  = X_test.loc[mask_tp, 'm_gg'].values
    mass_fn  = X_test.loc[mask_fn, 'm_gg'].values
    mass_bkg = X_test.loc[mask_bkg, 'm_gg'].values
    mass_sig = X_test.loc[mask_sig, 'm_gg'].values
    
    n_tp = len(mass_tp)
    n_fn = len(mass_fn)
    n_bkg = len(mass_bkg)
    n_sig = len(mass_sig)
    
    # --- Create the plot ---
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Histogram: Background (light gray, to show the full continuum)
    ax.hist(mass_bkg, bins=200, alpha=0.7,
            label=f'All Background (n={n_bkg:,})',
            color='#d3d3d3', edgecolor='gray', linewidth=0.5, histtype='stepfilled')
    
    # Histogram: Correctly identified signal (green, on top of everything)
    ax.hist(mass_sig, bins=200, alpha=0.3,
            label=f'Correct Signal (TP) (n={n_sig:,})',
            color='#2ca02c', edgecolor='darkgreen', linewidth=0.5, histtype='stepfilled')

    # Histogram: Wrongly classified signal (red, on top of background)
    #ax.hist(mass_fn, bins=200, alpha=0.7,
    #        label=f'Signal Misclassified (FN) (n={n_fn:,})',
    #        color='#d62728', edgecolor='darkred', linewidth=0.5, histtype='stepfilled')
    
    # Histogram: Correctly identified signal (green, on top of everything)
    #ax.hist(mass_tp, bins=200, alpha=0.8,
    #        label=f'Correct Signal (TP) (n={n_tp:,})',
    #        color='#2ca02c', edgecolor='darkgreen', linewidth=0.5, histtype='stepfilled')
    
    # Vertical line at π⁰ mass
    ax.axvline(x=135, color='black', linestyle='--', linewidth=1.5, label=r'$\pi^0$ mass (135 MeV)')
    
    # Labels and styling
    ax.set_xlabel(r'$m_{\gamma\gamma}$ [MeV/$c^2$]', fontsize=14)
    ax.set_ylabel('Events', fontsize=14)
    ax.set_title(plot_title or f'Mass Distribution Breakdown – {phys_ch}', fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

# =================================================================
# Plot BDT score signal breakdown (pair-wise)
# =================================================================
def plot_score_breakdown(y_test, y_pred, y_score, phys_ch="", plot_title=None):
    """
    Plot BDT score distribution with three categories:
        1. Correct Signal (True Positives)
        2. Signal Misclassified (False Negatives)
        3. All Background (True Negatives + False Positives)
    """
    print("Plotting score breakdown: TP, FN, All Background...")
    
    # --- Split the data ---
    mask_tp = (y_test == 1) & (y_pred == 1)
    mask_fn = (y_test == 1) & (y_pred == 0)
    mask_bkg = (y_test == 0)          # All background (TN + FP)
    mask_sig = (y_test == 1)          # All signal (TP + FN)
    
    score_tp = y_score[mask_tp]
    score_fn = y_score[mask_fn]
    score_bkg = y_score[mask_bkg]
    score_sig = y_score[mask_sig]
    
    n_tp = len(score_tp)
    n_fn = len(score_fn)
    n_bkg = len(score_bkg)
    n_sig = len(score_sig)
    
    # --- Create the plot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Histogram order: Background first, then misclassified, then correct
    #ax.hist([score_bkg, score_fn, score_tp], bins=100, alpha=0.7,
    #        label=[f'All Background (n={n_bkg:,})',
    #               f'Signal Misclassified (n={n_fn:,})',
    #               f'Correct Signal (n={n_tp:,})'],
    #        color=['#1f77b4', '#ff7f0e', '#2ca02c'],  # Blue, Orange, Green
    #        histtype='stepfilled', edgecolor='black', linewidth=0.5, stacked=False)
    
    # Histogram: Background
    ax.hist(score_bkg, bins=200, alpha=0.7,
            label=f'Background', # label=f'Background (n={n_bkg:,})'
            color='#d3d3d3', edgecolor='gray', linewidth=0.5, histtype='stepfilled')
    
    # Histogram: Signal
    ax.hist(score_sig, bins=200, alpha=0.3,
            label=f'Signal', # label=f'Signal (n={n_sig:,})'
            color='#2ca02c', edgecolor='darkgreen', linewidth=0.5, histtype='stepfilled')
    

    # Histogram: Wrongly classified signal (red, on top of background)
    #ax.hist(score_fn, bins=200, alpha=0.7,
    #        label=f'Signal Misclassified (FN) (n={n_fn:,})',
    #        color='#d62728', edgecolor='darkred', linewidth=0.5, histtype='stepfilled')
    
    # Histogram: Correctly identified signal (green, on top of everything)
    #ax.hist(score_tp, bins=200, alpha=0.8,
    #        label=f'Correct Signal (TP) (n={n_tp:,})',
    #        color='#2ca02c', edgecolor='darkgreen', linewidth=0.5, histtype='stepfilled')
    
    # Vertical line at typical cut threshold
    #ax.axvline(x=0.5, color='black', linestyle='--', linewidth=1.5, label='Cut threshold (0.5)')
    
    ax.set_xlabel('BDT Score', fontsize=16, fontweight='bold')
    ax.set_ylabel('Events', fontsize=16, fontweight='bold')
    #ax.set_title(plot_title, fontsize=16, fontweight='bold')
    #ax.set_title(plot_title or f'BDT Score Breakdown – {phys_ch}', fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')  # Log scale to see tails
    
    plt.tight_layout()
    return fig

# =================================================================
# Plot BDT score signal breakdown (event-wise)
# =================================================================
def plot_event_score_breakdown(event_df, score_col='max_proba', threshold=None, 
                               plot_title=None, phys_ch=""):
    """
    Plot event-level BDT score distribution:
        - Signal events (true_signal = 1)
        - Background events (true_signal = 0)
    Optionally mark a threshold.
    
    Args:
        event_df: DataFrame with columns 'true_signal' and score_col
        score_col: column name for the event score (e.g., 'max_proba')
        threshold: optional float to draw a vertical dashed line
        plot_title: custom title
        phys_ch: physics channel string (used if plot_title not given)
    """
    print("Plotting event-wise score breakdown...")
    
    sig_scores = event_df[event_df['true_signal'] == 1][score_col]
    bkg_scores = event_df[event_df['true_signal'] == 0][score_col]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Histogram: Background
    ax.hist(bkg_scores, bins=250, alpha=0.7,
            #label=f'Background events (n={len(bkg_scores):,})',
            label=f'Background events',
            color='#d3d3d3', edgecolor='gray', linewidth=0.5, histtype='stepfilled')
    
    # Histogram: Signal
    ax.hist(sig_scores, bins=250, alpha=0.5,
            #label=f'Signal events (n={len(sig_scores):,})',
            label=f'Signal events',
            color='#2ca02c', edgecolor='darkgreen', linewidth=0.5, histtype='stepfilled')
    
    # Optional threshold line
    if threshold is not None:
        ax.axvline(x=threshold, color='red', linestyle='--', linewidth=2,
                   label=f'Threshold = {threshold:.2f}')
    
    ax.set_xlim(-0.01, .8)
    #ax.set_ylim(-0.01, 1.01)
    ax.set_yscale('log')  # log scale to see tails

    ax.set_xlabel('BDT Score', fontsize=16, fontweight='bold')
    ax.set_ylabel('Events', fontsize=16, fontweight='bold')
    ax.set_title(plot_title , 
                 fontsize=16, fontweight='bold')
    #ax.set_title(plot_title or f'Event-wise Score Breakdown – {phys_ch}', 
    #             fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=16)   # <-- increased legend size
    ax.grid(True, alpha=0.3)
    
    
    plt.tight_layout()
    return fig

def plot_roc_improved(y_true, y_score, plot_title="", threshold=None):
    """
    Enhanced ROC curve with:
        - Larger fonts, thicker curve
        - Shaded area under the curve
        - Zoomed inset of the high-performance region
        - Metrics box
        - Optional operating point marker (if threshold is provided)
    """
    print("Plotting improved ROC curve...")
    
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    
    # --- Main figure ---
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    
    # Main ROC curve
    ax.plot(fpr, tpr, color='#1f77b4', lw=3, label=f'ROC AUC')
    
    # Shaded area under the curve
    ax.fill_between(fpr, tpr, alpha=0.15, color='#1f77b4')
    
    # Diagonal reference line
    ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Random classifier')
    
    # Axis labels and title
    ax.set_xlabel('False Positive Rate (FPR)', fontsize=18, fontweight='bold')
    ax.set_ylabel('True Positive Rate (TPR)', fontsize=18, fontweight='bold')
    #ax.set_title(plot_title or 'ROC Curve', fontsize=18, fontweight='bold')
    ax.set_title(plot_title, fontsize=18, fontweight='bold')   
    ax.legend(loc='upper right', fontsize=14, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Set axis limits
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)

    # ---------- MAIN AUC PLOT ----------
    ax.tick_params(axis='both', labelsize=18)   # or 16, 20, etc.
    
    # --- Metrics box ---
    idx_eff90 = np.argmin(np.abs(tpr - 0.90))
    fpr_eff90 = fpr[idx_eff90]
    bkg_rej_90 = 1 - fpr_eff90
    youden_idx = np.argmax(tpr - fpr)
    best_thresh = thresholds[youden_idx]
    best_tpr = tpr[youden_idx]
    best_fpr = fpr[youden_idx]
    
    textstr = (
        f"AUC: {roc_auc:.4f}\n"
        f"Bkg. rejection @ 90% S eff: {bkg_rej_90:.2%}\n"
        f"Optimal threshold (max TPR−FPR): {best_thresh:.3f}\n"
        f"    TPR = {best_tpr:.3f}, FPR = {best_fpr:.4f}"
    )
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
    ax.text(0.02, 0.02, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='bottom', horizontalalignment='left',
            bbox=props)
    
    # --- Operating point marking (if threshold is given) ---
    if threshold is not None:
        # Find the closest index to the given threshold
        idx = np.argmin(np.abs(thresholds - threshold))
        tpr_op = tpr[idx] #  coordinates of the operating point on the ROC curve
        fpr_op = fpr[idx]
        # Mark with a large red star
        ax.scatter(fpr_op, tpr_op, color='red', s=200, zorder=10, marker='*', edgecolors='black')
        # Annotate
        ax.text(fpr_op + 0.03, tpr_op - 0.1, 
                f'Cut = {threshold:.3f}\nTPR={tpr_op:.3f}\nFPR={fpr_op:.4f}',
                fontsize=10, color='red',
                bbox=dict(facecolor='white', edgecolor='red', boxstyle='round,pad=0.3'))
    
    # --- Inset zoom ---
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax, width="40%", height="40%", loc='lower right',
                       bbox_to_anchor=(0.05, 0.05, 0.9, 0.9),
                       bbox_transform=ax.transAxes)
    axins.set_xlim(-0.005, 0.2)
    axins.set_ylim(0.78, 1.01)
    axins.plot(fpr, tpr, color='#1f77b4', lw=2)
    axins.plot([0, 0.2], [0.8, 1], 'k--', lw=1, alpha=0.5)
    axins.grid(True, linestyle=':', alpha=0.4)
    axins.set_xlabel('FPR', fontsize=16)
    axins.set_ylabel('TPR', fontsize=16)
    from matplotlib.patches import Rectangle
    rect = Rectangle((0, 0.8), 0.2, 0.2, linewidth=1, edgecolor='gray',
                     facecolor='none', linestyle='--')
    ax.add_patch(rect)
    
    # If threshold is given, also show the point in the inset
    if threshold is not None:
        axins.scatter(fpr_op, tpr_op, color='red', s=80, zorder=10, marker='*')
    
    #plt.tight_layout()
    return fig

def plot_f1_comparison(strategies, thresholds, fpr_target=0.01, title=None):
    """
    Plot F1 score vs threshold for multiple strategies on the same axes.
    
    Parameters:
    -----------
    strategies : list of dicts, each with:
        - 'name' : str
        - 'f1_scores' : array-like
        - 'fpr_scores' : array-like (optional) – if provided, FPR is shown on secondary axis
        - 'color' : str (optional)
    thresholds : array-like (common for all)
    fpr_target : float, target FPR (default 0.01 = 1%)
    title : str, optional
    """
    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    
    # Primary y-axis: F1
    ax.set_xlabel('Threshold', fontsize=16, fontweight='bold')
    ax.set_ylabel('F1 Score', fontsize=16, fontweight='bold', color='black')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_xlim(0, 1.0)
    
    # Secondary axis for FPR (only if at least one strategy provides non-None fpr_scores)
    has_fpr = any('fpr_scores' in s and s['fpr_scores'] is not None for s in strategies)
    ax2 = None
    if has_fpr:
        ax2 = ax.twinx()
        ax2.set_ylabel('False Positive Rate (%)', fontsize=16, fontweight='bold', color='#d62728')
        ax2.tick_params(axis='y', labelcolor='#d62728')
    
    # Store handles for legend
    lines = []
    labels = []
    
    for i, strat in enumerate(strategies):
        name = strat['name']
        f1 = np.array(strat['f1_scores'])
        color = strat.get('color', f'C{i}')
        marker = strat.get('marker', 'o')
        
        # F1 curve
        line1, = ax.plot(thresholds, f1, marker=marker, linestyle='-', linewidth=2.5,
                         markersize=8, color=color, label=f'F1 ({name})')
        lines.append(line1)
        labels.append(f'F1 ({name})')
        
        # Mark maximum F1 point
        opt_idx = np.argmax(f1)
        opt_thr = thresholds[opt_idx]
        opt_f1 = f1[opt_idx]
        ax.plot(opt_thr, opt_f1, '*', color=color, markersize=18)
        ax.text(opt_thr + 0.02, opt_f1 - 0.1, f'{name}\n{opt_thr:.2f}', 
                fontsize=9, color=color)
        
        # Optional FPR curve – only if fpr_scores exists and is not None
        if 'fpr_scores' in strat and strat['fpr_scores'] is not None and ax2 is not None:
            fpr = np.array(strat['fpr_scores']) * 100
            line2, = ax2.plot(thresholds, fpr, marker='s', linestyle='--', linewidth=2,
                              markersize=7, color=color, alpha=0.7,
                              label=f'FPR ({name})')
            lines.append(line2)
            labels.append(f'FPR ({name})')
            
            # Mark threshold where FPR drops below target
            below_target = fpr < (fpr_target * 100)
            if np.any(below_target):
                idx_cross = np.argmax(below_target)
                thr_cross = thresholds[idx_cross]
                fpr_cross = fpr[idx_cross]
                ax2.plot(thr_cross, fpr_cross, 'g*', markersize=14)
                ax2.text(thr_cross + 0.02, fpr_cross + 1, f'FPR<1%', fontsize=9, color='green')
    
    # Add horizontal line for target FPR if secondary axis exists
    if ax2 is not None:
        ax2.axhline(y=fpr_target*100, color='green', linestyle='--', linewidth=1.5,
                    alpha=0.5, label=f'Target FPR = {fpr_target*100:.0f}%')
        ax2.tick_params(axis='both', labelsize=18)

    # Combine legends from both axes if needed
    if ax2 is not None:
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        all_lines = lines1 + lines2
        all_labels = labels1 + labels2
        # Remove duplicate label entries (keep only last occurrence)
        unique = {}
        for l, lab in zip(all_lines, all_labels):
            unique[lab] = l
        ax.legend(unique.values(), unique.keys(), loc='lower left', fontsize=12)
    else:
        ax.legend(loc='lower left', fontsize=16)
    
    plt.title(title, fontsize=18, fontweight='bold')
    
    #plt.title(title or 'F1 Score vs Threshold – Strategy Comparison', fontsize=18, fontweight='bold')
    return fig