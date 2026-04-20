import sys
import warnings

# Suppress pandas DBAPI2 warning
if not sys.warnoptions:
    warnings.simplefilter("ignore")

# Database layer (persistence)
import mysql.connector # mysql-connector-python
from mysql.connector import Error
import pandas as pd
from typing import List, Dict, Optional

class MySQLKLOEDB:
    #def __init__(self, host='localhost', database='kloe_bdt', 
    #             user='kloe_user', password='kloe_password'):
    def __init__(self, 
                 host='localhost', 
                 user='kloe_user', 
                 password='kloe_password'):
        # user = 'root', password = 'secret' or user = 'kloe_user', password = 'kloe_password'
        self.config = {
            'host': host,
            #'database': database,
            'user': user,
            'password': password,
            'raise_on_warnings': True
        }
        self.conn = None # Initialize attribute
        self.connect() # Immediately establish connection

    def connect(self):
        """
        Establish database connection
        """
        try:
            self.conn = mysql.connector.connect(**self.config)
            print("✅ Connected to MySQL")
        except Error as e:
            print(f"❌ Error: {e}")

    def create_database(self):
        """
        Create database if not exists
        """
        cursor = self.conn.cursor()
        try:
             # This creates a cursor object - a tool that lets you execute SQL commands and fetch results from the database.
            cursor.execute("CREATE DATABASE IF NOT EXISTS kloe_bdt")
            cursor.execute("USE kloe_bdt")
            print("Create database if not exists")
        except mysql.connector.errors.DatabaseError as e:
            if "atabase exists" in str(e) or "already exists" in str(e):
                # Database already exists
                print("ℹ️ Database already exists, continuing...")
                # Still need to USE the database
                cursor.execute("USE kloe_bdt")
            else:
                raise # Re-rasie if it's a different error
        
        self.conn.commit()

    def create_tables(self):
        """
        Create all required tables
        """
        cursor = self.conn.cursor()

        # Events table
        try:
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS events (
                    event_id INT AUTO_INCREMENT PRIMARY KEY,
                    run_number INT NOT NULL,
                    event_number INT NOT NULL,
                    is_signal BOOLEAN DEFAULT FALSE,
                    bdt_score DECIMAL(6,5),
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    INDEX idx_signal (is_signal),
                    INDEX idx_score (bdt_score)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
            """)
            print("\tEvents table is created!")
        except mysql.connector.errors.DatabaseError as e:
            if "already exists" in str(e):
                # Table already exists
                print("ℹ️ Table 'events' already exists, continuing...")
            else:
                raise e
        
        cursor.execute("""SHOW COLUMNS FROM events""") # Return rows and must fetch!
        columns = cursor.fetchall() # CONSUME the results

        #print columns
        for column in columns:
            print(f"\tColumn: {column[0]}, \tType: {column[1]}, \tNullable: {column[2]}")
        
        
        # Photon pairs table
        try:
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS photon_pairs (
                           pair_id INT AUTO_INCREMENT PRIMARY KEY,
                           event_id INT NOT NULL,
                           invariant_mass DECIMAL(10,3),
                           opening_angle DECIMAL(6,4),
                           energy_asymmetry DECIMAL(6,4),
                           energy_ratio DECIMAL(6,4),
                           energy_difference DECIMAL(10,3),
                           min_energy_angle DECIMAL(10,3),
                           asymmetry_angle DECIMAL(10,3),
                           bdt_prediction DECIMAL(6,5),
                           FOREIGN KEY (event_id) REFERENCES events(event_id) ON DELETE CASCADE,
                           INDEX idx_event (event_id),
                           INDEX idx_prediction (bdt_prediction)
                ) ENGINE=InnoDB 
            """)
            # ON DELETE CASCADE: If an event is deleted from the events table, all related photon pairs will be automatically deleted
            print("\tPhoton pair table is created!")
        except mysql.connector.errors.DatabaseError as e:
            if "already exists" in str(e):
                # Table already exists
                print("ℹ️ Table 'photon pairs' already exists, continuing...")
            else:
                raise e
            
        cursor.execute("""SHOW COLUMNS FROM photon_pairs""") # Return rows 
        columns = cursor.fetchall() 
        for column in columns:
            print(f"\tColumn: {column[0]}, \tType: {column[1]}, \tNullable: {column[2]}")

        # Features table (for ML traning)

        # Model metadata table
        
        self.conn.commit()
        print(f"✅ Tables created")

    def insert_event(self, run_number: int, event_number: int, bdt_score: float, is_signal: bool) -> int:
        """
        Insert event and return event_id
        """
        print("Insert event and return event_id")
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT INTO events (run_number, 
                       event_number, 
                       bdt_score, 
                       is_signal)
            VALUES (%s, %s, %s, %s)               
        """, (run_number, event_number, bdt_score, is_signal)) # %s is a placeholder
        self.conn.commit()
        return cursor.lastrowid
    
    def insert_photon_pair(self, event_id: int, features: Dict):
        """
        Insert photon pair features
        """
        print("Insert photon pair features")
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT INTO photon_pairs (
                       event_id, 
                       invariant_mass,
                       opening_angle,
                       energy_asymmetry,
                       energy_ratio,
                       energy_difference,
                       min_energy_angle,
                       asymmetry_angle,
                       bdt_prediction
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
        """, (
            event_id,
            features.get('invariant_mass'),  
            features.get('opening_angle'),     
            features.get('energy_asymmetry'),
            features.get('energy_ratio'),      
            features.get('energy_difference'),
            features.get('min_energy_angle'),
            features.get('asymmetry_angle'),
            features.get('bdt_prediction')
        ))
        self.conn.commit()

    def query_signal_events(self, min_score: float = 0.8, limit: int = 100) -> pd.DataFrame:
        """
        Query high-confidence signal events
        """
        query = """
            SELECT e.event_id,
                e.run_number, 
                e.event_number, 
                e.bdt_score,
                e.is_signal,
                p.invariant_mass,
                p.opening_angle,
                p.energy_asymmetry,
                p.energy_ratio,
                p.energy_difference,
                p.min_energy_angle,
                p.asymmetry_angle,
                p.bdt_prediction
            FROM events e
            JOIN photon_pairs p ON e.event_id = p.event_id
            WHERE e.bdt_score > %s AND e.is_signal = 1
            ORDER BY e.bdt_score DESC
            LIMIT %s
        """
        return pd.read_sql(query, self.conn, params=[min_score, limit])

    def close(self):
        if self.conn: # Is there a connection object?
            self.conn.close()

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

