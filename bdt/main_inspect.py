import joblib
import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)
print(project_root)

from plots import plot_compr_hist, plot_var, plot_feature_pairs, plot_feature_target
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from config import (
    DATA_DIR, PLOTS_INSPECT_DIR
)

# Inspect features and correlations

if __name__ == '__main__':

    #DATA_DIR = 'dataset' 
    #input_data_dir = os.path.join(project_root, f'analysis/{DATA_DIR}')
    input_data_dir = DATA_DIR
    phys_map = joblib.load(os.path.join(input_data_dir, f'phys_map.pkl'))
    print(phys_map)

    # Create output folder
    #plot_dir = os.path.join(project_root, 'analysis/plots')
    plot_dir = PLOTS_INSPECT_DIR
    
    # With this (always fresh) plot_dir
    import shutil
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir, exist_ok=True)

    ## Check if true pi0 events add up for all channels in phys_map
    nb_evnt_sum = 0
    nb_evnt_combined = 0
    for data_type, info in phys_map.items():
        #br_nm = data_type.split(';')[0]
        br_title = info['br_title']
        category = info['category']

        file_path = os.path.join(input_data_dir, f'all_df_{data_type}.pkl')
        if not os.path.exists(file_path):
            print(f"❌ File not found: {file_path}")
            print(f"   This channel was likely skipped due to low statistics (<10 events)")
            continue
        
        df_tmp = joblib.load(os.path.join(input_data_dir, f'all_df_{data_type}.pkl'))
        nb_evnt_tmp = len(df_tmp)
        #print(f"key: {data_type}, category: {category}")

        if category == 'combined':
            nb_evnt_combined = nb_evnt_tmp
            print(f"{data_type}: {nb_evnt_tmp}")
        else:
            nb_evnt_sum += nb_evnt_tmp
            print(f"{data_type}: {nb_evnt_tmp}")

    print(f"nb evnt combined: {nb_evnt_combined}, sum: {nb_evnt_sum}")

    #============================================================
    # INSPECTING DATA. ONE AT THE TIME
    #============================================================
    
    #phys_ch = ['TISR3PI_SIG', 'signal']
    #phys_ch = ['TOMEGAPI', 'background']
    #phys_ch = ['TKPM', 'background']
    #phys_ch = ['TKSL', 'background']
    #phys_ch = ['TRHOPI', 'background']
    #phys_ch = ['TBKGREST', 'background']
    #phys_ch = ['TEEG', 'background']

    phys_ch = ['TCOMB', 'combined']
    
    for data_type, info in phys_map.items():
        br_nm = data_type.split(';')[0]
        br_title = info['br_title']
        category = info['category']
        print(f"Inspecting dataset {br_nm}; {br_title}; {category}")  
        #print(info)

        if (br_nm  == phys_ch[0]and category == phys_ch[1]):
            print(f'{info}')

            ## Check if dataset exists
            
            
            try:
                all_df = joblib.load(os.path.join(input_data_dir, f'all_df_{data_type}.pkl'))
                print("Loaded existing pickle file")
            except FileNotFoundError:
                print("Pickle file not found, quit")
                sys.exit(1) 

            ## Load dataset
            pos_df = all_df[all_df['is_signal'] == 1]
            neg_df = all_df[all_df['is_signal'] == 0]   
            #print(all_df.columns)
            #print(pos_df.head(5))
            #print(neg_df.head(5))

            pi0_all_df = joblib.load(os.path.join(input_data_dir, f'pi0_all_df_{data_type}.pkl'))
            pi0_pos_df = pi0_all_df[pi0_all_df['is_pi0'] == 1]
            pi0_neg_df = pi0_all_df[pi0_all_df['is_pi0'] == 0]
            #print(pi0_all_df.columns)
            pi0_mass = pi0_pos_df['m_gg'].tolist()

            ## Inpect kinematic variables
            kine_features = ['Br_betapi0', 'Br_ppIM', 'Br_angle_pi0gam12', 'Br_deltaE', 'Br_m3pi', 'Br_lagvalue_min_7C', 'is_signal']
            kine_all_df = all_df[kine_features]
            kine_pos_df = kine_all_df[kine_all_df['is_signal'] == 1]
            kine_neg_df = kine_all_df[~(kine_all_df['is_signal'] == 1)]    
            #
            #print(kine_all_df['Br_deltaE'].describe())
            #print(kine_all_df.columns)  


            ## * Plot gamma 4-momentum 
            df_set = [all_df, pos_df, neg_df]
            drop_columns = ['event', 'Br_betapi0', 'Br_ppIM', 'Br_angle_pi0gam12', 'Br_deltaE', 'Br_m3pi', 
                            'Br_lagvalue_min_7C', 'Br_recon_indx', 'Br_bkg_indx', 'is_signal', 'true_pi0_pair']
            
            fig_compr_hist = plot_compr_hist(df_set, drop_columns,
                                             3, 100, 
                                             rf"$\gamma$ 4-momentum ({br_title})") # Photon 4-momentum comparison plot
            fig_compr_hist.savefig(f'{plot_dir}/Photon_4-momentum_compr_{br_nm}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_compr_hist)
            #print(all_df.columns.tolist())

            ## * Plot pi0 features
            drop_columns = ['event', 'pair_id', 'is_pi0']
            pi0_df_set = [pi0_all_df, pi0_pos_df, pi0_neg_df]
            
            fig_compr_hist = plot_compr_hist(pi0_df_set, drop_columns,
                                             2, 100, 
                                             rf"$\pi^{0}$ Candidates ({br_title})") # Pi0 comparison plot
            fig_compr_hist.savefig(f'{plot_dir}/Pi0_compr_{br_nm}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_compr_hist)  

            ## * Plot kine. var of selection cuts: chi2, gg-opening-angle, deltaE, ppIM, betapi0
            drop_columns = ['is_signal']
            kine_df_set = [kine_all_df, kine_pos_df, kine_neg_df]
            fig_compr_hist = plot_compr_hist(kine_df_set, drop_columns,
                                             2, 100, 
                                             rf"Selection Cuts ({br_title})") # Kine. var. comparison plot
            fig_compr_hist.savefig(f'{plot_dir}/Kine_compr_{br_nm}.png', dpi=300, bbox_inches='tight')
            plt.close(fig_compr_hist)
            #print(kine_all_df.head(5)) 

            
            if category == 'signal' or category == 'combined': # Only plots for signal channels
                print(f"!!!!!!!!!!!!!!!!!!!!!!!!!! {br_nm}, {category}")
                ## * Plot true pi0 mass
                fig_var = plot_var(pi0_mass, 'm_gg', br_title)
                fig_var.savefig(f'{plot_dir}/pos_pi0_mass_{br_nm}.png', dpi=300, bbox_inches='tight')
                plt.close(fig_var)

                ## * Plot kine. var of selection cuts correlations
                drop_columns = ['']
                fig_feature_pairs = plot_feature_pairs(kine_all_df, drop_columns,
                                                       rf"Selection Cuts Correlations (Signal=Blue, Background=Red) ({br_title})", 
                                                       "is_signal")
                fig_feature_pairs.savefig(f'{plot_dir}/SC_correlation_{br_nm}.png', dpi=300, bbox_inches='tight')
                plt.close(fig_feature_pairs.fig)  

                ## * Plot pi0 feature-feature correlations
                print(pi0_all_df.head(5))
                drop_columns = ['event', 'cos_theta', 'pair_id']
                fig_feature_pairs = plot_feature_pairs(pi0_all_df, drop_columns,
                                                       rf"$\pi^{0}$ Candidates Feature-feature (Signal=Blue, Background=Red) ({br_title})", 
                                                       "is_pi0")
                fig_feature_pairs.savefig(f'{plot_dir}/FF_correlation_{br_nm}.png', dpi=300, bbox_inches='tight')
                plt.close(fig_feature_pairs.fig)

                ## * Plot pi0 feature-target correlations
                features = ['m_gg', 'opening_angle', 'cos_theta', 'E_asym', 'e_min_x_angle', 
                            'E1', 'E2', 'E3', 'asym_x_angle', 'E_diff', 'is_pi0']
                target_corr = pi0_all_df[features].corr()['is_pi0'].drop('is_pi0') #.sort_values(ascending=False)
                #sorted_by_abs = target_corr.abs().sort_values(ascending=False)
                fig_feature_target = plot_feature_target(target_corr, rf'Feature Importance: Correlation with true $\pi^{0}$ ({br_title})')
                fig_feature_target.savefig(f'{plot_dir}/feature_target_correlation_{br_nm}.png', dpi=300, bbox_inches='tight')
                plt.close(fig_feature_target)
    
