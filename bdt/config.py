"""
Simple configuration management for pi0-gamma-pairing project
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import json
# Project root
#PROJECT_ROOT = Path(__file__).parent.absolute()
PROJECT_ROOT = (Path(__file__).parent.absolute()).parent.absolute()



# Data directories
DATA_DIR = PROJECT_ROOT / 'dataset_bdt'
PLOTS_INSPECT_DIR = PROJECT_ROOT / 'plots_inspect'

#MODEL_DIR = PROJECT_ROOT / 'training' / 'models'
#PLOT_DIR_VAL = PROJECT_ROOT / 'validation' / 'plots'
#PLOT_DIR_APP = PROJECT_ROOT / 'application' / 'plots'



# Initializing configuration
REQUIRED_BR = ['Br_E1', 'Br_px1', 'Br_py1', 'Br_pz1', 
                        'Br_E2', 'Br_px2', 'Br_py2', 'Br_pz2', 
                        'Br_E3', 'Br_px3', 'Br_py3', 'Br_pz3',
                        'Br_m3pi', 'Br_lagvalue_min_7C', 'Br_deltaE',
                        'Br_angle_pi0gam12', 'Br_ppIM', 'Br_betapi0',
                        'Br_recon_indx', 'Br_bkg_indx']

# =================================================================
# CREATE DATASET
# =================================================================

def create_dataset(df, category): # For photon 4-momentum
    print(f'\n✅ Creating dataset ...')
    # Define photon 4-momentum
    br_nm = ['Br_E1', 'Br_px1', 'Br_py1', 'Br_pz1', 
             'Br_E2', 'Br_px2', 'Br_py2', 'Br_pz2', 
             'Br_E3', 'Br_px3', 'Br_py3', 'Br_pz3',
             'Br_m3pi', 'Br_lagvalue_min_7C', 'Br_deltaE',
             'Br_angle_pi0gam12', 'Br_ppIM', 'Br_betapi0',
             'Br_recon_indx', 'Br_bkg_indx']
    
    # ADDED: Check which columns actually exist
    available_cols = [col for col in br_nm if col in df.columns]
    if len(available_cols) < len(br_nm):
        print(f"  Note: Using {len(available_cols)}/{len(br_nm)} available columns")
    
    # Selection cut, ensure physical region
    cut_region = (df['Br_lagvalue_min_7C'] < 100) if 'Br_lagvalue_min_7C' in df.columns else pd.Series(True, index=df.index)
    phys_region = ((df['Br_betapi0'] < 1) & (df['Br_betapi0'] > 0)) if 'Br_betapi0' in df.columns else pd.Series(True, index=df.index)

    #df = df[(df['Br_lagvalue_min_7C'] < 100) & (df['Br_betapi0'] < 1) & (df['Br_betapi0'] > 0)][br_nm]
    # Apply filters and create a proper copy
    df_filtered = df[cut_region & phys_region][available_cols].copy()
    df = df_filtered  # Now df is a clean copy

    # Create all_df, pos_df, neg_df for signal and background events
    if len(available_cols): # Check para length and br_nm length are consistent
        
        if category == 'signal':
            print(f"Creating all_df for {category} {df.columns}...")

            # ADDED: Check if classification columns exist
            if 'Br_recon_indx' in df.columns and 'Br_bkg_indx' in df.columns:
                pos_df = df[(df['Br_recon_indx'] == 2) & (df['Br_bkg_indx'] == 1)][available_cols]
                neg_df = df[~((df['Br_recon_indx'] == 2) & (df['Br_bkg_indx'] == 1))][available_cols]
            else:
                print(f"  Warning: Missing classification columns, treating all as background")
                pos_df = pd.DataFrame(columns=available_cols)
                neg_df = df[available_cols]

            # True positive
            nb_pos = [i for i in range(len(pos_df))]  
            if len(pos_df) > 0:
                pos_df = pos_df.copy()
                pos_df.insert(0, 'event', nb_pos)  # Add event column
                pos_df['is_signal'] = 1 # Add is_signal column
                pos_df['true_pi0_pair'] = [(0, 1)] * len(nb_pos)  # Add true_pi0_pair column
            
            # True negative
            nb_neg = [i for i in range(len(neg_df))]
            if len(neg_df) > 0:
                neg_df = neg_df.copy()
                neg_df.insert(0, 'event', nb_neg) # Add event column 
                neg_df['is_signal'] = 0 # Add is_signal column
                neg_df['true_pi0_pair'] = [(-1, -1)] * len(nb_neg) # Add true_pi0_pair column

            # Combine pos + neg dataset and shuffling
            dfs_to_concat = []
            if len(pos_df) > 0:
                dfs_to_concat.append(pos_df)
            if len(neg_df) > 0:
                dfs_to_concat.append(neg_df)
            
            if dfs_to_concat:
                all_df = pd.concat(dfs_to_concat, ignore_index=True)
                all_df = all_df.sample(frac=1).reset_index(drop=True)
            else:
                all_df = pd.DataFrame()
      
        elif category == 'background':
            print(f"Creating pho4mom_all_df for {category} {df.columns} ...")

            all_df = df.copy()
            nb_all_df = [i for i in range(len(all_df))]  
            all_df.insert(0, 'event', nb_all_df)  # Add event column
            all_df['is_signal'] = 0 # Add is_signal column
            all_df['true_pi0_pair'] = [(-1, -1)] * len(all_df) # Add true_pi0_pair column
        else: # combined or others
            raise ValueError("Only sig and bkg allow!")
    else: 
        print("No dataset other than signal or bkg or combined is expected!")
        raise ValueError("Array length mismatch - cannot proceed")

    # pi0 features for ML learning
    # MODIFIED: Only create pairs if we have data
    if len(all_df) > 0:
        pi0_all_df = prepare_3photon_paris(all_df)
    else:
        pi0_all_df = pd.DataFrame()
        
    return all_df, pi0_all_df

# =================================================================
# DATA SPLITTING
# =================================================================

def data_splitting(all_df):
    print(f'\n✅ Splitting dataset ...')

    # First split: separate test set 20%
    all_df_trainval, all_df_test = train_test_split(
        all_df, test_size=0.2, random_state=42
    )

    # Second split: separte validation from training (20% of total)
    all_df_train, all_df_val = train_test_split(
        all_df_trainval, test_size=0.25, random_state=42
    )

    print(f"Training events:   {len(all_df_train)} ({len(all_df_train)/len(all_df)*100:.1f}%)")
    print(f"Validation events: {len(all_df_val)} ({len(all_df_val)/len(all_df)*100:.1f}%)")
    print(f"Test events:       {len(all_df_test)} ({len(all_df_test)/len(all_df)*100:.1f}%)")

    #print("CREATING PAIRS FROM EACH SPLIT")

    pair_train = prepare_3photon_paris(all_df_train)
    pair_val = prepare_3photon_paris(all_df_val)
    pair_test = prepare_3photon_paris(all_df_test)

    # Verify event column preservation
    for name, pair_df in [('Train', pair_train), ('Val', pair_val), ('Test', pair_test)]:
        if 'event' not in pair_df.columns:
            print(f"❌ CRITICAL: 'event' column missing in {name} pairs!")
            print(f"   Columns: {pair_df.columns.tolist()}")
        else:
            n_events = pair_df['event'].nunique()
            n_pairs = len(pair_df)
            print(f"✅ {name}: {n_pairs} pairs from {n_events} events ({n_pairs/n_events:.2f} pairs/event)")

    #print(f"Training pairs:   {len(pair_train)}")
    #print(f"Validation pairs: {len(pair_val)}")
    #print(f"Test pairs:       {len(pair_test)}")

    ## Step 3: Prepare features for training
    features = ['m_gg', 'opening_angle', 'cos_theta', 'E_asym', 'e_min_x_angle', 'E1', 'E2', 'E3', 'asym_x_angle', 'E_diff']
    
    # ADDED: Check which features actually exist
    available_features = [f for f in features if f in pair_train.columns]
    if len(available_features) < len(features):
        print(f"  Warning: Using {len(available_features)}/{len(features)} available features")
    
    X_train = pair_train[available_features] if len(pair_train) > 0 else pd.DataFrame()
    y_train = pair_train['is_pi0'] if 'is_pi0' in pair_train.columns else pd.Series()
    X_val = pair_val[available_features] if len(pair_val) > 0 else pd.DataFrame()
    y_val = pair_val['is_pi0'] if 'is_pi0' in pair_val.columns else pd.Series()
    X_test = pair_test[available_features] if len(pair_test) > 0 else pd.DataFrame()
    y_test = pair_test['is_pi0'] if 'is_pi0' in pair_test.columns else pd.Series()

    print(f"\n✅ Data ready for training:")
    print(f"  X_train shape: {X_train.shape}, y_train shape: {y_train.shape}")
    print(f"  X_val shape:   {X_val.shape}, y_val shape:   {y_val.shape}")
    print(f"  X_test shape:  {X_test.shape}, y_test shape:  {y_test.shape}")

    return all_df_train, all_df_val, all_df_test, X_train, y_train, X_val, y_val, X_test, y_test, pair_train, pair_val, pair_test

# =================================================================
# Prpare 3 photon paris
# =================================================================
def prepare_3photon_paris(df_events):
    """
    Convert 3-photon events into training paris with EXACT pi0 invariant masses.

    Assumes your DataFrame has columns:
    E1, px1, py1, pz1,  # OR E1, pt1, eta1, phi1
    E2, px2, py2, pz2,
    E3, px3, py3, pz3,
    is_signal, true_pi0_pair
    """

    #print("="*50)
    #print("Creating photon pairs with EXACT quantities: invariant masses ...")
    #print("="*50)

    pairs = []

    for _, evt in df_events.iterrows():
        # Get 4-vector for all 3 photons
        # ADAPT THIS TO YOUR EXACT COLUMN NAMES
        photons = [# [E, px, py, pz]
            np.array([evt.Br_E1, evt.Br_px1, evt.Br_py1, evt.Br_pz1]), 
            np.array([evt.Br_E2, evt.Br_px2, evt.Br_py2, evt.Br_pz2]),
            np.array([evt.Br_E3, evt.Br_px3, evt.Br_py3, evt.Br_pz3])
        ]

        # All 3 possible pairs.
        #pair_indices = [(0,1), (0,2), (1,2)]
        pair_indices = [(0,1), (2,0), (1,2)]

        #print(f"{pair_indices},{type(pair_indices)}, {type(df_events)}")
        #print(photon)

        for i, j in pair_indices:

            # Extract energy values
            e1 = photons[i][0]
            e2 = photons[j][0]

            # Check for missing values in individual floats
            e1_missing = pd.isna(e1) or (isinstance(e1, float) and np.isnan(e1))
            e2_missing = pd.isna(e2) or (isinstance(e2, float) and np.isnan(e2))

            # Handle missing values in calculations
            energy_threshold = 5  # 5 MeV energy cut, aviod small e1+e2 causes E_asym close to 1
            if e1_missing or e2_missing or e1 < energy_threshold or e2 < energy_threshold:
                # If either energy is missing, set derived values to NaN
                e_asym = np.nan
                e_ratio = np.nan
                e_diff = np.nan
                # Skip angle calculations that depend on energy
                theta = np.nan
                cos_theta = np.nan
                e_min_x_angle = np.nan
                asym_x_angle = np.nan
                mass = np.nan
            else:
                # Calculate EXACT invariant mass from 4-vectors
                mass = inv_mass_4vector(photons[i], photons[j])
                #print(f"mass = {mass}")

                #print(f"paired indicds ({i}, {j}), unpaired index {unpaired_idx}")

                # Opening angle
                p1_mag = np.sqrt(np.maximum(0., photons[i][1]**2 + photons[i][2]**2 + photons[i][3]**2))
                p2_mag = np.sqrt(np.maximum(0., photons[j][1]**2 + photons[j][2]**2 + photons[j][3]**2))
                dot_product = photons[i][1] * photons[j][1] + photons[i][2] * photons[j][2] + photons[i][3] * photons[j][3]
                cos_theta = dot_product / (p1_mag * p2_mag + 1e-10) # 1e-10 avoid divide zero
                cos_theta = np.clip(cos_theta, -1, 1) # Forced to physical range
                theta = np.arccos(np.clip(cos_theta, -1, 1)) # Range (0, pi)
                #print(f"p1_mag = {p1_mag}, p2_mag = {p2_mag}")
                #print(f"dot_product = {dot_product}, cos_theta = {cos_theta}")
                #print(f"{np.clip(cos_theta, -1, -1)}")
                #print(f"theta = {theta}")

                # Energy asymmetry
                e_asym = np.abs(e1 - e2) / (e1 + e2 + 1e-10)
                e_asym = np.clip(e_asym, 0, 1)  # Force to physical range
                #print(f"p_asym = {e_asym}")

                # Energy ratio
                e_ratio = min(e1, e2) / max(e1, e2)

                # Energy diff.
                e_diff = np.abs(e1 - e2)
                
                # Minimum energy angle
                e_min_x_angle = min(e1, e2) * theta

                # Energy asymmetry angle
                asym_x_angle = e_asym * theta

            # Find the unpaired photon index (the one not in {i,j})
            unpaired_idx = [k for k in range(3) if k not in [i, j]][0]
            e3 = photons[unpaired_idx][0]
            e3_missing = pd.isna(e3) or (isinstance(e3, float) and np.isnan(e3))

            # Is this the correct pi0 pair? (require truth info)
            is_pi0 = 0
            if 'is_signal' in evt and evt.is_signal == 1:
                if 'true_pi0_pair' in evt:
                    is_pi0 = 1 if (i, j) == evt.true_pi0_pair else 0
                else:
                    # If you don't have exact pair truth,
                    # assume the pair with mass closest to 0.135 GeV is correct
                    is_pi0 = 1 if abs(mass - 0.135) < 0.015 else 0

            # Data insertion
            pairs.append({
                'event': evt.event,
                'pair_id': f"{evt.event}_{i}{j}",
                'm_gg': mass,
                'opening_angle': theta, # Opening angle in radians
                'cos_theta': cos_theta,
                'E_asym': e_asym,
                'e_min_x_angle': e_min_x_angle,
                'E1': e1,
                'E2': e2,
                'E3': e3,
                'asym_x_angle': asym_x_angle,
                'E_diff': e_diff,
                #'E_ratio': e_ratio,
                'is_pi0': is_pi0
            })
        
    pi0_all_df = pd.DataFrame(pairs)

    # Deal with extreame values: E_asym, asym_x_angle
    #
    E_asym_values = [pair['E_asym'] for pair in pairs]
    E_asym_df = pd.DataFrame(E_asym_values, columns=['E_asym'])
    #print(E_asym_df.describe())

    #
    asym_x_angle_values = [pair['asym_x_angle'] for pair in pairs]
    asym_x_angle_df = pd.DataFrame(asym_x_angle_values, columns=['asym_x_angle'])
    #print(asym_x_angle_df.describe())
    
        
    #print('m_gg', pi0_all_df['m_gg'].describe())
    #print(pi0_all_df['E_asym'].describe())

    return pi0_all_df

# =================================================================
# Invariant mass 4-vector
# =================================================================
def inv_mass_4vector(p1, p2):
    """
    Calculate diphoton invariant mass from two photon 4-vectors.

    Args:
        p1, p2: Arrays/lists of [E, px, py, pz] or [E, pt, eta, phi]

    Returns:
        Invariant mass
    """

    if len(p1) == 4 and len(p2) == 4:
        e = p1[0] + p2[0]
        px = p1[1] + p2[1]
        py = p1[2] + p2[2]
        pz = p1[3] + p2[3]
        mass_squared = e**2 - (px**2 + py**2 + pz**2)
        # Ensure non-negative before sqrt
        return np.sqrt(np.maximum(0., mass_squared))
        #if (mass_squared < 0):
            #print(f"mass_squared    {mass_squared}")
        #return np.sqrt(mass_squared)
    else:   
        # Use your experiment's Lorentz vector class
        # (e.g., ROOT.TLorentzVector, vector.obj, etc.)
        return (p1 + p2).M()

# =================================================================
# PATCH FUNCTION
# =================================================================
def patched_get_basescore(model):
    config_str = model.get_booster().save_config()
    config = json.loads(config_str)
    base_score_str = config["learner"]["learner_model_param"]["base_score"]
    # Remove brackets if present
    base_score_str = base_score_str.strip('[]')
    return float(base_score_str)