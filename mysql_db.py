# Database layer (persistence)
import mysql.connector # mysql-connector-python
from mysql.connector import Error
import pandas as pd
from typing import List, Dict, Optional

class MySQLKLOEDB:
    def __init__(self, host='localhost', database='kloe_bdt', 
                 user='kloe_user', password='kloe_password'):
        # user = 'root', password = 'secret' or user = 'kloe_user', password = 'kloe_password'
        self.config = {
            'host': host,
            'database': database,
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
                cursor = self.conn.cursor()
                cursor.execute("USE kloe_bdt")
            else:
                raise # Re-rasie if it's a different error
        
        self.conn.commit()

    def create_tables(self):
        """
        Create all required tables
        """
        cursor = self.conn.cursor()

        try:
            # Events table
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
            print(f"\tColumn: {column[0]}, \tType: {column[1]}, \tNullable: {columns[2]}")
        
        self.conn.commit()
          
    def insert_event(self, run_number: int, event_number: int, bdt_score: float, is_signal: bool) -> int:
        """
        Insert event and return event_id
        """
        print("Insert event and return event_id")
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT INTO events (run_number, event_number, bdt_score, is_signal)
            VALUES (%s, %s, %s, %s)               
        """, (run_number, event_number, bdt_score, is_signal)) # %s is a placeholder
        self.conn.commit()
        return cursor.lastrowid

    def close(self):
        if self.conn: # Is there a connection object?
            self.conn.close()

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

# Usage
with MySQLKLOEDB() as db:
    db.create_database()
    db.create_tables()