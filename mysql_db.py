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
    def __init__(self, host='localhost', database='kloe_bdt', 
                 user='kloe_user', password='kloe_password'):
        self.host = host
        self.database = database
        self.user = user
        self.password = password
        self.conn = None
        
        # Try to connect with database first
        if not self._try_connect_with_db():
            # If that fails, connect without database
            self._connect_without_db()
    
    def _try_connect_with_db(self):
        """Try to connect with the database specified"""
        try:
            self.conn = mysql.connector.connect(
                host=self.host,
                database=self.database,
                user=self.user,
                password=self.password,
                raise_on_warnings=True
            )
            print(f"✅ Connected to MySQL database: {self.database}")
            return True
        except Error as e:
            if "Unknown database" in str(e):
                print(f"ℹ️ Database '{self.database}' doesn't exist yet")
                return False
            else:
                print(f"❌ Connection Error: {e}")
                return False
    
    def _connect_without_db(self):
        """Connect without specifying a database"""
        try:
            self.conn = mysql.connector.connect(
                host=self.host,
                user=self.user,
                password=self.password,
                raise_on_warnings=True
            )
            print(f"✅ Connected to MySQL server (no database selected)")
        except Error as e:
            print(f"❌ Connection Error: {e}")
            raise
    
    def connect(self):
        """Establish database connection (legacy method for compatibility)"""
        if self.conn is None or not self.conn.is_connected():
            self._connect_without_db()
    
    def create_database(self):
        """Create database if not exists"""
        if self.conn is None:
            self.connect()
        
        cursor = self.conn.cursor()
        try:
            cursor.execute(f"CREATE DATABASE IF NOT EXISTS {self.database}")
            cursor.execute(f"USE {self.database}")
            self.conn.commit()
            print(f"✅ Database '{self.database}' created/selected")
        except Error as e:
            print(f"❌ Database creation error: {e}")
            raise
        finally:
            cursor.close()

    def create_tables(self):
        """Create all required tables"""
        if self.conn is None:
            raise Exception("No database connection")
        
        cursor = self.conn.cursor()
        try:
            cursor.execute(f"USE {self.database}")
            
            # Drop and recreate events table with correct precision
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS events (
                    event_id INT AUTO_INCREMENT PRIMARY KEY,
                    run_number INT NOT NULL,
                    event_number INT NOT NULL,
                    is_signal BOOLEAN DEFAULT FALSE,
                    bdt_score FLOAT,  
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    INDEX idx_signal (is_signal),
                    INDEX idx_score (bdt_score)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
            """)
            print("✅ Events table ready")
            
            # Photon pairs table
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
                    bdt_prediction FLOAT,  # Changed from DECIMAL(6,5) to DECIMAL(7,6)
                    FOREIGN KEY (event_id) REFERENCES events(event_id) ON DELETE CASCADE,
                    INDEX idx_event (event_id),
                    INDEX idx_prediction (bdt_prediction)
                ) ENGINE=InnoDB
            """)
            print("✅ Photon pairs table ready")
            
            self.conn.commit()
            print("✅ All tables created/verified")
            
        except Error as e:
            print(f"❌ Table creation error: {e}")
            raise
        finally:
            cursor.close()

    def insert_event(self, run_number: int, event_number: int, bdt_score: float, is_signal: bool) -> int:
        """Insert event and return event_id"""
        if self.conn is None:
            raise Exception("No database connection")
        
        cursor = self.conn.cursor()
        try:
            cursor.execute("""
                INSERT INTO events (run_number, event_number, bdt_score, is_signal)
                VALUES (%s, %s, %s, %s)               
            """, (run_number, event_number, bdt_score, is_signal))
            self.conn.commit()
            event_id = cursor.lastrowid
            return event_id
        except Error as e:
            print(f"❌ Error inserting event: {e}")
            raise
        finally:
            cursor.close()
    
    def insert_photon_pair(self, event_id: int, features: Dict):
        """Insert photon pair features"""
        if self.conn is None:
            raise Exception("No database connection")
        
        cursor = self.conn.cursor()
        try:
            cursor.execute("""
                INSERT INTO photon_pairs (
                    event_id, invariant_mass, opening_angle, energy_asymmetry,
                    energy_ratio, energy_difference, min_energy_angle, asymmetry_angle, bdt_prediction
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
        except Error as e:
            print(f"❌ Error inserting photon pair: {e}")
            raise
        finally:
            cursor.close()

    def query_signal_events(self, min_score: float = 0.8, limit: int = 100) -> pd.DataFrame:
        """Query high-confidence signal events"""
        query = """
            SELECT e.event_id, e.run_number, e.event_number, e.bdt_score,
                   e.is_signal, p.invariant_mass, p.opening_angle, p.energy_asymmetry,
                   p.energy_ratio, p.energy_difference, p.min_energy_angle,
                   p.asymmetry_angle, p.bdt_prediction
            FROM events e
            JOIN photon_pairs p ON e.event_id = p.event_id
            WHERE e.bdt_score > %s AND e.is_signal = 1
            ORDER BY e.bdt_score DESC 
            LIMIT %s
        """
        return pd.read_sql(query, self.conn, params=[min_score, limit])

    def get_training_data(self, limit: int = 10000) -> pd.DataFrame:
        """Export data for model retraining"""
        query = """
            SELECT 
                p.invariant_mass, p.opening_angle, p.energy_asymmetry,
                p.energy_ratio, p.energy_difference, p.min_energy_angle,
                p.asymmetry_angle, e.is_signal as label
            FROM photon_pairs p
            JOIN events e ON p.event_id = e.event_id
            WHERE e.is_signal IS NOT NULL
            LIMIT %s
        """
        return pd.read_sql(query, self.conn, params=[limit])

    def close(self):
        if self.conn and self.conn.is_connected():
            self.conn.close()
            print("✅ Database connection closed")

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()