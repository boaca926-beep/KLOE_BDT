-- init_db.sql
-- this script will be automatically run when the MySQL container activates for the first time

CREATE DATABASE IF NOT EXISTS kloe_bdt;
USE kloe_bdt;

CREATE TABLE IF NOT EXISTS model_metadata (
    id INT AUTO_INCREMENT PRIMARY KEY,
    model_name VARCHAR(100),
    model_version VARCHAR(50),
    is_active BOOLEAN DEFAULT TRUE,
    training_date DATETIME,
    n_features INT,
    feature_names TEXT,
    auc_score FLOAT,
    accuracy FLOAT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS events (
    event_id INT AUTO_INCREMENT PRIMARY KEY,
    run_number INT,
    event_number INT,
    bdt_score FLOAT,
    is_signal BOOLEAN,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS photon_pairs (
    pair_id INT AUTO_INCREMENT PRIMARY KEY,
    event_id INT,
    invariant_mass FLOAT,
    opening_angle FLOAT,
    energy_asymmetry FLOAT,
    energy_ratio FLOAT,
    energy_difference FLOAT,
    min_energy_angle FLOAT,
    asymmetry_angle FLOAT,
    bdt_prediction FLOAT,
    FOREIGN KEY (event_id) REFERENCES events(event_id) ON DELETE CASCADE
);