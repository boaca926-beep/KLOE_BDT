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
        try:
            cursor = self.conn.cursor() # This creates a cursor object - a tool that lets you execute SQL commands and fetch results from the database.
            cursor.execute("CREATE DATABASE IF NOT EXISTS kloe_bdt")
            cursor.execute("USE kloe_bdt")
            self.conn.commit()
            print("Create database if not exists")
        except mysql.connector.errors.DatabaseError as e:
            if "" in str(e):
                # Database already exists
                print("ℹ️ Database already exists, continuing...")
                # Still need to USE the database
                ursor = self.conn.cursor()
                cursor.execute("USE kloe_bdt")
            else:
                raise # Re-rasie if it's a different error

    def create_tables(self):
        """
        Create all required tables
        """
        cursor = self.conn.cursor()

        # Events table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS events (
                event_id INT AUTO_INCREMENT PRIMARY KEY,
                run_number INT NOT NULL,
                event_number INT NOT NULL,
                is_signal BOOLEAN DEFAULT FALSE,
                bdt_score DECIMAL(6,5),
                create_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                INDEX idx_signal (is_signal),
                INDEX idx_score (bdt_score)
            ) ENGINE=InnoDB DEFAULT CHARSET = utf8mb4
        """)
        print("\tEvents table is created!")

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