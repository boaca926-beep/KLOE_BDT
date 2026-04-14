# Training configs
import pandas as pd
import numpy as np
import random

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