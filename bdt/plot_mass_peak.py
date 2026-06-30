import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from prediction import find_best_pi0_candidate
from config import DATA_DIR, MODEL_DIR
from scipy.optimize import curve_fit

# ========== CREATE OUTPUT DIRECTORY ==========
plot_dir = "../plot_mass_peak"
import shutil
if os.path.exists(plot_dir):
    print(f"🧹 Clearing existing directory: {plot_dir}")
    shutil.rmtree(plot_dir)
os.makedirs(plot_dir, exist_ok=True)
print(f"✅ Plots directory created: {plot_dir}")

# ========== LOAD TEST SET ==========
data_type = 'TCOMB'
input_data_dir = DATA_DIR

all_df_test = joblib.load(os.path.join(input_data_dir, f'all_df_test_{data_type}.pkl'))
X_test = joblib.load(os.path.join(input_data_dir, f'X_test_{data_type}.pkl'))
y_test = joblib.load(os.path.join(input_data_dir, f'y_test_{data_type}.pkl'))
model = joblib.load(os.path.join(MODEL_DIR, f'pi0_classifier_model_{data_type}.pkl'))

print(f"Loaded {len(all_df_test)} photons, {len(X_test)} pairs")

# ========== COMPUTE PER-PHOTON MASS AND SCORE ==========
results = []
for idx, evt in all_df_test.iterrows():
    photons = [
        np.array([evt.Br_E1, evt.Br_px1, evt.Br_py1, evt.Br_pz1]),
        np.array([evt.Br_E2, evt.Br_px2, evt.Br_py2, evt.Br_pz2]),
        np.array([evt.Br_E3, evt.Br_px3, evt.Br_py3, evt.Br_pz3]),
    ]
    best_pair, score, mass = find_best_pi0_candidate(photons, model)
    results.append({
        'mass': mass,
        'score': score,
        'true_signal': evt['is_signal'],
    })

df = pd.DataFrame(results)
print(f"Computed mass and score for {len(df)} photons")

# ========== SELECT OPTIMAL THRESHOLD ==========
threshold = 0.15
selected = df[df['score'] > threshold]
print(f"Selected {len(selected)} photons (threshold = {threshold})")

# ========== PLOT 1: BDT SCORE DISTRIBUTION ==========
plt.figure(figsize=(8,6))
sig_scores = df[df['true_signal']==1]['score']
bkg_scores = df[df['true_signal']==0]['score']
plt.hist(sig_scores, bins=50, alpha=0.5, label='Signal', density=True)
plt.hist(bkg_scores, bins=50, alpha=0.5, label='Background', density=True)
plt.xlabel('BDT Score (mean strategy)')
plt.ylabel('Density')
plt.legend()
plt.title(r'BDT Score Distribution for Test Set')
plt.savefig(os.path.join(plot_dir, 'bdt_score_distribution.png'), dpi=300)
print(f"Saved: {os.path.join(plot_dir, 'bdt_score_distribution.png')}")

# ========== PLOT 2: MASS COMPARISON (BEFORE vs AFTER BDT CUT) ==========
plt.figure(figsize=(8,6))
plt.hist(df['mass'], bins=100, range=(100, 180), alpha=0.5, label='Before BDT cut (All)', density=True, histtype='step')
plt.hist(selected['mass'], bins=100, range=(100, 180), alpha=0.5, label=f'After BDT cut (score > {threshold})', density=True, histtype='step')
plt.xlabel(r'$m_{\gamma\gamma}$ (MeV/$c^2$)')
plt.ylabel('Density')
plt.legend()
plt.title(r'$\pi^0$ Mass Peak Before and After BDT Selection')
plt.savefig(os.path.join(plot_dir, 'pi0_mass_comparison.png'), dpi=300)
print(f"Saved: {os.path.join(plot_dir, 'pi0_mass_comparison.png')}")

# ========== PLOT 3: MASS DISTRIBUTION WITH SEPARATE SIGNAL/BACKGROUND ==========
# Use total selected mass for fitting
mass_data = selected['mass'].values
fit_range = (110, 160)
bins = np.linspace(fit_range[0], fit_range[1], 160)
hist, bin_edges = np.histogram(mass_data, bins=bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Double Gaussian fit on total
def double_gaussian(x, amp1, mean1, sigma1, amp2, mean2, sigma2):
    return amp1 * np.exp(-0.5 * ((x - mean1) / sigma1)**2) + \
           amp2 * np.exp(-0.5 * ((x - mean2) / sigma2)**2)

p0 = [400, 135, 5, 100, 135, 10]
bounds = ([0, 100, 1, 0, 100, 1], [10000, 150, 15, 5000, 150, 20])

try:
    popt, pcov = curve_fit(double_gaussian, bin_centers, hist, p0=p0, bounds=bounds, maxfev=5000)
    amp1, mean1, sigma1, amp2, mean2, sigma2 = popt
    # reduced chi2
    errors = np.sqrt(hist + 1e-9)
    residuals = hist - double_gaussian(bin_centers, *popt)
    chi2_val = np.sum((residuals / errors)**2)
    ndof = len(bin_centers) - len(popt)
    chi2_red = chi2_val / ndof
    # signal yield
    x_range = np.linspace(fit_range[0], fit_range[1], 1000)
    dx = (fit_range[1] - fit_range[0]) / 1000
    signal_hist = double_gaussian(x_range, *popt)
    signal_yield = np.sum(signal_hist) * dx

    print("\n=== Double Gaussian Fit (Total Selected Mass) ===")
    print(f"Main Gaussian: amp = {amp1:.1f}, mean = {mean1:.3f} MeV, sigma = {sigma1:.3f} MeV")
    print(f"Second Gaussian: amp = {amp2:.1f}, mean = {mean2:.3f} MeV, sigma = {sigma2:.3f} MeV")
    print(f"Reduced chi² = {chi2_red:.3f} (ndof={ndof})")
    print(f"Signal yield (integral) = {signal_yield:.1f} events")
except Exception as e:
    print("Fit failed:", e)
    popt = None

# Separate signal and background after cut
sig_selected = selected[selected['true_signal'] == 1]['mass']
bkg_selected = selected[selected['true_signal'] == 0]['mass']

fig, ax = plt.subplots(figsize=(10, 7))
# Total histogram (step)
ax.hist(mass_data, bins=bins, range=fit_range, histtype='step', color='black', linewidth=1.5, label='Total (selected)')
# Signal histogram (filled, blue, transparent)
ax.hist(sig_selected, bins=bins, range=fit_range, alpha=0.5, color='blue', label='True Signal')
# Background histogram (filled, red, transparent)
ax.hist(bkg_selected, bins=bins, range=fit_range, alpha=0.5, color='red', label='True Background')

if popt is not None:
    x_plot = np.linspace(fit_range[0], fit_range[1], 500)
    ax.plot(x_plot, double_gaussian(x_plot, *popt), 'r-', linewidth=2, label='Double Gaussian fit (total)')
    
    # Format fit parameters into a compact, two‑row text box
    # Row 1: main and tail parameters
    row1 = (r'Main: $\mu_1={:.2f}$ MeV, $\sigma_1={:.2f}$ MeV').format(mean1, sigma1, mean2, sigma2)
    row2 = (r'Tail: $\mu_2={:.2f}$ MeV, $\sigma_2={:.2f}$ MeV').format(mean1, sigma1, mean2, sigma2)
    row3 = r'$\chi^2/\mathrm{{ndof}} = {:.3f}$'.format(chi2_red)
    textstr = row1 + '\n' + row2 + '\n' + row3

    # Place text in top‑right corner (0.95, 0.95) with right‑alignment
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

ax.set_xlabel(r'$m_{\gamma\gamma}$ (MeV/$c^2$)')
ax.set_ylabel('Entries / bin')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_title(r'$\pi^0$ Mass Peak: Total, Signal, and Background after BDT Cut')
plt.savefig(os.path.join(plot_dir, 'pi0_mass_separated.png'), dpi=300)
print(f"Saved: {os.path.join(plot_dir, 'pi0_mass_separated.png')}")

print("\n✅ All plots generated in", plot_dir)