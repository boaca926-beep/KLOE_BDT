import numpy as np

# =================================================================
# For a NEW 3-photon event, pick the best pi0 candidate
# =================================================================
def find_best_pi0_candidate(photon_4vectors, model):
    """
    photon_4vectors: list of 3 arrays, each [E, px, py, pz] or [E, pt, eta, phi]
    Returns: (best_pair_indices, probability, mass)
    """
    pairs = [(0,1), (0,2), (1,2)]
    candidates = []

    for i, j in pairs:
        # Calculate EXACT quantities from 4-vectors
        mass = inv_mass_4vector(photon_4vectors[i], photon_4vectors[j])
        #print(f"mass    {mass}; photon1_4vectors    {photon_4vectors[i]}; photon2_4vectors  {photon_4vectors[j]}")

        # True opening angle
        p1 = photon_4vectors[i]
        p2 = photon_4vectors[j]
        p1_mag = np.sqrt(np.maximum(0., p1[1]**2 + p1[2]**2 + p1[3]**2))
        p2_mag = np.sqrt(np.maximum(0., p2[1]**2 + p2[2]**2 + p2[3]**2))
        dot_product = p1[1] * p2[1] + p1[2] * p2[2] + p1[3] * p2[3]
        cos_theta = dot_product / (p1_mag * p2_mag + 1e-10) # 1e-10 avoid divide zero
        theta = np.arccos(np.clip(cos_theta, -1, 1))
        #print(f"theta:  {theta}")

        # Energy asymmetry
        e1, e2 = p1[0], p2[0]
        e_asym = np.abs(e1 - e2) / (e1 + e2 + 1e-10)
        # print(f"e1: {e1}; e2: {e2}; asym: {asym}")

        # Minimum energy angle
        e_min_x_angle = min(e1, e2) * theta

        # Find the unpaired photon index (the one not in {i,j})
        unpaired_idx = [k for k in range(3) if k not in [i, j]][0]
        e3 = photon_4vectors[unpaired_idx][0]

        # Energy asymmetry angle
        asym_x_angle = e_asym * theta

        # Energy diff.
        e_diff = np.abs(e1 - e2)

        features_values = [mass, theta, cos_theta, e_asym, e_min_x_angle, e1, e2, e3, asym_x_angle, e_diff]

        # Predict
        proba = model.predict_proba([features_values])[0, 1]

        candidates.append({
            'pair': (i, j),
            'score': proba,
            'mass': mass,
            'theta': theta
        })

        # Return the best candidate
        best = max(candidates, key=lambda x: x['score']) # best: The entire dictionary with the highest score
        # candidates: A list of dictionaries, each with a 'score' key
        # max(): Built-in Python function that finds the maximum value
        # key=lambda x: x['score']: Tells max() to use the 'score' value for comparison
        return best['pair'], best['score'], best['mass']
    
# =================================================================
# Invariant mass of 2 gamma 4-vector
# =================================================================
def inv_mass_4vector(p1, p2):
    """
    Calculate diphoton invariant mass from two photon 4-vectors.

    Args:
        p1, p2: Arrays/lists of [E, px, py, pz] or [E, pt, eta, phi]

    Returns:
        Invariant mass in GeV
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