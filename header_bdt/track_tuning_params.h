#ifndef TRACK_TUNING_PARAMS_H
#define TRACK_TUNING_PARAMS_H

// ============================================================
// POST-FIT TRACK SMEARING TUNING PARAMETERS
// ============================================================
// extra_smear: constant additional smearing added to tracks
// after the kinematic fit (to match data RMS)
//
// Tuning procedure:
// 1. Start with extra_smear = 0.0
// 2. Run on a small sub-sample (1000-5000 events)
// 3. Measure ppIM RMS for each extra_smear value
// 4. Interpolate to find value that gives RMS = 70 MeV
// 5. Update this value and re-run full analysis
// ============================================================

const double EXTRA_SMEAR = 0.010;   // <--- TUNE THIS VALUE

#endif
