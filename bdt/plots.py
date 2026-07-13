# Plotting functions
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, auc, accuracy_score
from sklearn.metrics import confusion_matrix, classification_report

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
        print("⚠️ Warning: No feature columns to plot!")
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
    #plt.show(block=False)

    # Get detailed metrics
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))

    return fig