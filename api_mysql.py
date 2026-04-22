# api_mysql.py - Updated to handle 10 features and fix data truncation

from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel
from typing import List, Optional, Dict
import joblib
import pandas as pd
import numpy as np
from mysql_db import MySQLKLOEDB
import os

app = FastAPI(title="KLOE BDT API (MySQL)")

# Load model
try:
    model = joblib.load("models/pi0_classifier_model_TCOMB.pkl")
    print(f"✅ Model loaded! Expected features: {model.n_features_in_}")
    if hasattr(model, 'feature_names_in_'):
        print(f"   Features: {model.feature_names_in_}")
except Exception as e:
    print(f"❌ Failed to load model: {e}")
    model = None

# Database dependency
def get_db():
    db = MySQLKLOEDB(
        host=os.getenv('MYSQL_HOST', 'localhost'),  
        user=os.getenv('MYSQL_USER', 'kloe_user'),  
        password=os.getenv('MYSQL_PASSWORD', 'kloe_password')  
    )
    try:
        cursor = db.conn.cursor()
        cursor.execute("USE kloe_bdt")
        yield db
    except Exception as e:
        print(f"Database error: {e}")
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    finally:
        db.close()

# Define the 10 features your model expects
class PhotonPairFeatures(BaseModel):
    m_gg: float
    opening_angle: float
    cos_theta: float
    E_asym: float
    e_min_x_angle: float
    E1: float
    E2: float
    E3: float
    asym_x_angle: float
    E_diff: float

class PredictionRequest(BaseModel):
    run_number: int
    event_number: int
    photon_pairs: List[PhotonPairFeatures]

@app.post("/predict-and-save")
async def predict_and_save(
    request: PredictionRequest,
    db: MySQLKLOEDB = Depends(get_db)
):
    """Predict and immediately save to MySQL"""
    if model is None:
        raise HTTPException(status_code=500, detail="Model not loaded")
    
    if not request.photon_pairs:
        raise HTTPException(status_code=400, detail="No photon pairs provided")
    
    results = []

    for pair in request.photon_pairs:
        # Convert to DataFrame with correct feature order
        feature_dict = pair.dict()
        df = pd.DataFrame([feature_dict])
        
        # Ensure columns are in the order the model expects
        if hasattr(model, 'feature_names_in_'):
            df = df[model.feature_names_in_]
        
        # Make prediction
        probabilities = model.predict_proba(df)
        score = float(probabilities[:, 1][0])
        is_signal = score > 0.35

        # Save to database
        event_id = db.insert_event(
            request.run_number,
            request.event_number,
            score,
            is_signal
        )
        
        # Calculate and sanitize values for database
        # Energy ratio (E1/E2) - ensure it's between 0-1 and rounded to 4 decimal places
        energy_ratio = pair.E1 / pair.E2 if pair.E2 > 0 else 0
        energy_ratio = max(0, min(1, energy_ratio))  # Clamp between 0-1
        energy_ratio = round(energy_ratio, 4)  # Round to 4 decimal places
        
        # Other values that need rounding
        invariant_mass = round(pair.m_gg, 3)
        opening_angle = round(pair.opening_angle, 4)
        energy_asymmetry = round(pair.E_asym, 4)
        energy_difference = round(pair.E_diff, 3)
        min_energy_angle = round(pair.e_min_x_angle, 3)
        asymmetry_angle = round(pair.asym_x_angle, 4)
        
        # Prepare features for database (7 columns)
        photon_features = {
            'invariant_mass': invariant_mass,
            'opening_angle': opening_angle,
            'energy_asymmetry': energy_asymmetry,
            'energy_ratio': energy_ratio,
            'energy_difference': energy_difference,
            'min_energy_angle': min_energy_angle,
            'asymmetry_angle': asymmetry_angle,
            'bdt_prediction': round(score, 6)
        }
        
        db.insert_photon_pair(event_id, photon_features)

        results.append({
            'event_id': event_id,
            'bdt_score': score,
            'is_signal': is_signal,
            'energy_ratio': energy_ratio  # For debugging
        })

    return {
        "status": "saved", 
        "run_number": request.run_number,
        "event_number": request.event_number,
        "total_pairs": len(request.photon_pairs),
        "successful": len(results),
        "results": results
    }

@app.get("/events/signal")
async def get_signal_events(
    min_score: float = 0.8,
    limit: int = 100,
    db: MySQLKLOEDB = Depends(get_db)
):
    """Retrieve high-quality signal events"""
    try:
        df = db.query_signal_events(min_score, limit)
        
        if df.empty:
            return {"message": "No signal events found", "events": []}
        
        return df.to_dict(orient='records')
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Query failed: {str(e)}")

@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "api": "KLOE BDT API with MySQL",
        "version": "2.0",
        "status": "running",
        "model_loaded": model is not None,
        "n_features": model.n_features_in_ if model else 0,
        "endpoints": {
            "POST /predict-and-save": "Make predictions and save to database",
            "GET /events/signal": "Get high-confidence signal events",
            "GET /health": "Health check",
            "GET /docs": "Swagger documentation"
        }
    }

@app.get("/health")
async def health_check(db: MySQLKLOEDB = Depends(get_db)):
    """Health check endpoint"""
    try:
        cursor = db.conn.cursor()
        cursor.execute("SELECT DATABASE()")
        current_db = cursor.fetchone()[0]
        
        return {
            "status": "healthy",
            "database": f"connected to {current_db}",
            "model": "loaded" if model else "not loaded",
            "features_expected": model.n_features_in_ if model else 0
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Health check failed: {str(e)}")