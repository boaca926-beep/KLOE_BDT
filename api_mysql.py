# API layer (business logic)
from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel
from typing import List, Optional
import joblib
import pandas as pd
#from mysql_db import MySQLKLOEDB

"""
Application Programming Interface (API)
"""

app = FastAPI(title="KLOE BDT API (MySQL)")

# Load model
model = joblib.load("models/pi0_classifier_model_TCOMB.pkl")
print("Model is loaded!")