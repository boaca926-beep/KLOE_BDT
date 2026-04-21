## 📊 MySQL Implementation

## 📚 Table of Contents
- [Overview](#-overview)

## 💡 Overview
```text
KLOE_BDT Deployment/
├── PROJECT_SWD.md        # ⭐ Most important
├── mysql_db.py           # Database layer (persistence)
├── api_mysql.py          # API layer (business logic)
├── initi_bdt.py          # Initialize 'kloe_bdt' database, and fill in a few test events
├── models/
│   └── pi0_classifier_model_TCOMB.pkl  # Trained ML model
├── requirements.txt          # Dependencies
├── Dockerfile                # Dockerfile
├── docker-compose.mysql.yml  # Container orchestration
└── PRJOECT_SWD.md            # Documentation
```
<div align="center">
<img src="plots_ref/deployment_schema.png" width="500" alt="ROC & AUC"/>
<br/>
<em></em>

</div>

## 🚀 Work Flow

### 1. Install and Setup

```bash
# Ubuntu/Debian
sudo apt install mysql-server mysql-client
sudo mysql_secure_installation

# Docker (easiest)
docker run --name mysql-kloe -e MYSQL_ROOT_PASSWORD=secret -p 3306:3306 -d mysql:8
```


### 2. Python MySQL Connector

```python
# mysql_db.py
import mysql.connector
from mysql.connector import Error
import pandas as pd
from typing import List, Dict, Optional

class MySQLKLOEDB:
    def __init__(self, host='localhost', database='kloe_bdt', user='root', password='secret'):
        self.config = {
            'host': host,
            'database': database,
            'user': user,
            'password': password,
            'raise_on_warnings': True
        }
        self.conn = None
        self.connect()
    
    def connect(self):
        """Establish database connection"""
        try:
            self.conn = mysql.connector.connect(**self.config)
            print("✅ Connected to MySQL")
        except Error as e:
            print(f"❌ Error: {e}")
    
    def create_database(self):
        """Create database if not exists"""
        cursor = self.conn.cursor()
        cursor.execute("CREATE DATABASE IF NOT EXISTS kloe_bdt")
        cursor.execute("USE kloe_bdt")
        self.conn.commit()
    
    def create_tables(self):
        """Create all required tables"""
        cursor = self.conn.cursor()
        
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
                bdt_prediction DECIMAL(6,5),
                FOREIGN KEY (event_id) REFERENCES events(event_id) ON DELETE CASCADE,
                INDEX idx_event (event_id),
                INDEX idx_prediction (bdt_prediction)
            ) ENGINE=InnoDB
        """)
        
        # Features table (for ML training)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS features (
                feature_id INT AUTO_INCREMENT PRIMARY KEY,
                event_id INT NOT NULL,
                feature_name VARCHAR(50),
                feature_value DECIMAL(10,4),
                FOREIGN KEY (event_id) REFERENCES events(event_id) ON DELETE CASCADE,
                INDEX idx_event_feature (event_id, feature_name)
            ) ENGINE=InnoDB
        """)
        
        # Model metadata table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS model_metadata (
                id INT AUTO_INCREMENT PRIMARY KEY,
                model_name VARCHAR(100),
                model_version VARCHAR(20),
                training_date TIMESTAMP,
                auc_score DECIMAL(6,5),
                accuracy DECIMAL(6,5),
                threshold DECIMAL(6,5)
            )
        """)
        
        self.conn.commit()
        print("✅ Tables created")
    
    def insert_event(self, run_number: int, event_number: int, bdt_score: float, is_signal: bool) -> int:
        """Insert event and return event_id"""
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT INTO events (run_number, event_number, bdt_score, is_signal)
            VALUES (%s, %s, %s, %s)
        """, (run_number, event_number, bdt_score, is_signal))
        self.conn.commit()
        return cursor.lastrowid
    
    def insert_photon_pair(self, event_id: int, features: Dict):
        """Insert photon pair features"""
        cursor = self.conn.cursor()
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
    
    def query_signal_events(self, min_score: float = 0.8, limit: int = 100) -> pd.DataFrame:
        """Query high-confidence signal events"""
        query = """
            SELECT e.event_id, e.run_number, e.event_number, e.bdt_score,
                   p.invariant_mass, p.opening_angle, p.energy_asymmetry
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
        if self.conn:
            self.conn.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

```

### 3. Test MySQL Connector
```python
from mysql_db import MySQLKLOEDB
import mysql.connector # mysql-connector-python
from mysql.connector import Error
import pandas as pd


def check_database_exists():
    """Check if kloe_bdt database exists"""

    try:
        # Connect without specifying database
        conn = mysql.connector.connect(
            host='localhost',
            user='kloe_user',
            password='kloe_password'
        )
        cursor = conn.cursor()

        # Get all databases
        cursor.execute("SHOW DATABASES")
        databases = cursor.fetchall()

        # Check if kloe_bdt exists
        for db in databases:
            if db[0] == 'kloe_bdt':
                print(f"Database {db[0]} EXISTS")

                # To see tables in the database
                cursor.execute(f"SHOW TABLES FROM {db[0]}")
                tables = cursor.fetchall()

                if tables:
                    print(f"Tables in {db[0]}")
                    for table in tables:
                        table_name = table[0]
                        print(f"\t- {table_name}; Table content:")
                        # Show table structrue
                        cursor.execute(f"SELECT * FROM {db[0]}.{table_name}")
                        content = cursor.fetchall()
                        if content:
                            # Get column names
                            cursor.execute(f"DESCRIBE {db[0]}.{table_name}")
                            columns = cursor.fetchall()
                            col_names = [col[0] for col in columns]

                            # Print column headers
                            print(f"\t\t Columns: {', '.join(col_names)}, {type(columns)}, columns[0]: {type(columns[0])}")

                            # Print each row
                            for i, row in enumerate(content, 1):
                                print(f"\t\t Row {i}: {row}")


                cursor.close()
                conn.close()
                return True
        
        print(f"Database 'kloe_bdt' does NOT exist")
        cursor.close()
        conn.close()
        return False
    
    except Error as e:
        print(f"Error checking database: {e}")
        return False
    
def drop_database_via_root():
    """
    Drop database using root (more privileges)
    """
    try:
        conn = mysql.connector.connect(
            host='localhost',
            user='root',
            password='secret' # root password
        )
        cursor = conn.cursor()
        cursor.execute('DROP DATABASE IF EXISTS kloe_bdt')
        conn.commit()
        print("Database 'kloe_bdt' dropped via root!")
        cursor.close()
        conn.close()
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False
    

if __name__ == '__main__':

    kloe_db_status = check_database_exists()
    #print(f"{kloe_db_status}")

    if kloe_db_status:
        print("\n All data in kloe_bdt will be deleted!")
        confirm = input("Type 'DROP' to delete: ")

        if confirm == 'DROP':
            print("Dropping 'kloe_bdt' database")

            drop_database_via_root()

        else:

            print("Operation cancelled")
        
        check_database_exists()
    
    else:
        #print("Database 'kloe_bdt' doesn't exist!")
        confirm = input("Type 'YES' to initialize 'kloe_bdt' database : ")

        if confirm == "YES":

            # initialize database
            with MySQLKLOEDB() as db:
                db.create_database()
                db.create_tables()

                # Insert data
                r''' Single events with photon_pair
                event_id = db.insert_event(12345, 67890, 0.95, True)
                db.insert_photon_pair(event_id,{
                    'invariant_mass': 135.2,
                    'opening_angle': 0.85,
                    'energy_asymmetry': 0.12,
                    'energy_ratio': 0.88,
                    'energy_difference': 15.3,
                    'min_energy_angle': 42.1,
                    'asymmetry_angle': 0.10,
                    'bdt_prediction': 0.96
                })
                '''

                # Insert multiple events with their photon pairs
                events_with_pairs = [
                    {
                        'event': (12345, 67890, 0.95, True),
                        'photon_pair': {
                            'invariant_mass': 135.2,
                            'opening_angle': 0.85,
                            'energy_asymmetry': 0.12,
                            'energy_ratio': 0.88,
                            'energy_difference': 15.3,
                            'min_energy_angle': 42.1,
                            'asymmetry_angle': 0.10,
                            'bdt_prediction': 0.96
                        }
                    },
                    {
                        'event': (12345, 67891, 0.92, True),
                        'photon_pair': {
                            'invariant_mass': 134.8,
                            'opening_angle': 0.87,
                            'energy_asymmetry': 0.08,
                            'energy_ratio': 0.92,
                            'energy_difference': 10.2,
                            'min_energy_angle': 38.5,
                            'asymmetry_angle': 0.07,
                            'bdt_prediction': 0.94
                        }
                    },
                    {
                        'event': (12345, 67892, 0.88, True),
                        'photon_pair': {
                            'invariant_mass': 135.5,
                            'opening_angle': 0.82,
                            'energy_asymmetry': 0.15,
                            'energy_ratio': 0.85,
                            'energy_difference': 18.7,
                            'min_energy_angle': 45.2,
                            'asymmetry_angle': 0.13,
                            'bdt_prediction': 0.89
                        }
                    }
                ]
                
                for data in events_with_pairs:
                    event_id = db.insert_event(*data['event']) # Unpack tuple *, ** for dirctionaries
                    db.insert_photon_pair(event_id, data['photon_pair'])
                    print(f"Inserted event {event_id} with BDT score {data['event'][2]}")

                db_status = check_database_exists()
                #print(db_status)

                # Query results
                # See all events
                bdt_score=0.9

                print(f"SELECTED EVENTS WITH BDT_SCORE > {bdt_score}")
                all_events = pd.read_sql("SELECT * FROM events", db.conn)
                print("\t=== All Events ===")
                print(all_events)

                # See all photon pairs
                all_pairs = pd.read_sql("SELECT * FROM photon_pairs", db.conn)
                print("\t=== All Photon Pairs ===")
                print(all_pairs)

                # Join query
                
                results = db.query_signal_events(min_score=0.9)
                print(results.to_string())

        else:
            print("❌ Quit")
```


### 4. FastAPI + MySQL Integration
```python
# api_mysql.py
from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel
from typing import List, Optional, Dict
import joblib
import pandas as pd
from mysql_db import MySQLKLOEDB

app = FastAPI(title="KLOE BDT API (MySQL)")

# Load model
model = joblib.load("models/pi0_classifier_model_TCOMB.pkl")

# Database dependency
def get_db():
    db = MySQLKLOEDB(password='secret')
    try:
        yield db
    finally:
        db.close()

class PredictionRequest(BaseModel):
    run_number: int
    event_number: int
    photon_pairs: List[Dict]

@app.post("/predict-and-save")
async def predict_and_save(
    request: PredictionRequest,
    db: MySQLKLOEDB = Depends(get_db)
):
    """Predict and immediately save to MySQL"""
    results = []
    
    for pair in request.photon_pairs:
        # Make prediction
        df = pd.DataFrame([pair])
        score = float(model.predict_proba(df)[:, 1][0])
        is_signal = score > 0.35
        
        # Save to database
        event_id = db.insert_event(
            request.run_number, 
            request.event_number, 
            score, 
            is_signal
        )
        db.insert_photon_pair(event_id, {**pair, 'bdt_prediction': score})
        
        results.append({
            'event_id': event_id,
            'bdt_score': score,
            'is_signal': is_signal
        })
    
    return {"status": "saved", "results": results}

@app.get("/events/signal")
async def get_signal_events(
    min_score: float = 0.8,
    limit: int = 100,
    db: MySQLKLOEDB = Depends(get_db)
):
    """Retrieve high-quality signal events"""
    df = db.query_signal_events(min_score, limit)
    return df.to_dict(orient='records')
```

### 5. Docker Compose with MySQL
```yaml
# docker-compose.mysql.yml
version: '3.8'
services:
  mysql:
    image: mysql:8
    container_name: kloe-mysql
    environment:
      MYSQL_ROOT_PASSWORD: secret
      MYSQL_DATABASE: kloe_bdt
      MYSQL_USER: ml_user
      MYSQL_PASSWORD: ml_password
    ports:
      - "3306:3306"
    volumes:
      - mysql_data:/var/lib/mysql
    command: --default-authentication-plugin=mysql_native_password
    healthcheck:
      test: ["CMD", "mysqladmin", "ping", "-h", "localhost"]
      timeout: 10s
      retries: 5

  phpmyadmin:
    image: phpmyadmin/phpmyadmin
    container_name: kloe-phpmyadmin
    environment:
      PMA_HOST: mysql
      PMA_PORT: 3306
    ports:
      - "8080:80"
    depends_on:
      - mysql

  api:
    build: .
    container_name: kloe-api
    ports:
      - "8000:8000"
    depends_on:
      mysql:
        condition: service_healthy
    environment:
      MYSQL_HOST: mysql
      MYSQL_USER: ml_user
      MYSQL_PASSWORD: ml_password
      MYSQL_DATABASE: kloe_bdt
    command: uvicorn api_mysql:app --host 0.0.0.0 --port 8000 --reload

volumes:
  mysql_data:
```

### 6. MySQL-Specific Optimizations for Your Use Case
```sql
-- Indexes for faster queries (add to create_tables)
CREATE INDEX idx_bdt_score ON events(bdt_score);
CREATE INDEX idx_signal_score ON events(is_signal, bdt_score);

-- Partition by run_number for large datasets
ALTER TABLE events PARTITION BY RANGE (run_number) (
    PARTITION p2020 VALUES LESS THAN (202100),
    PARTITION p2021 VALUES LESS THAN (202200),
    PARTITION p2022 VALUES LESS THAN (202300),
    PARTITION p_future VALUES LESS THAN MAXVALUE
);

-- Stored procedure for batch insertion
DELIMITER $$
CREATE PROCEDURE BatchInsertEvents(
    IN run_num INT,
    IN start_event INT,
    IN end_event INT
)
BEGIN
    DECLARE i INT DEFAULT start_event;
    WHILE i <= end_event DO
        INSERT INTO events (run_number, event_number, created_at)
        VALUES (run_num, i, NOW());
        SET i = i + 1;
    END WHILE;
END$$
DELIMITER ;
```

### 7. MySQL on Your CV
```markdown
**Technical Skills**
- Databases: MySQL (advanced), SQLite, query optimization, indexing
- API Development: FastAPI, RESTful design, async endpoints
- ML Engineering: XGBoost, GPU acceleration, model deployment
- DevOps: Docker, docker-compose, MySQL replication
```
## Note
### MySQL 
#### Installation
**Check if MySQL is Installed**
```bash
# Check package status
dpkg -l | grep mysql-server
# Or
apt list --installed | grep mysql-server

# Check binary location
which mysql
# Output: /usr/bin/mysql (if installed)

# Check version
mysql --version
# Output: mysql  Ver 8.0.36-0ubuntu0.22.04.1 for Linux on x86_64

# Check all MySQL packages
dpkg -l | grep -E "mysql|mariadb"
```

**Check if MySQL is Running (Activated)**
```bash
# Systemctl (most common)
sudo systemctl status mysql
# Look for: active (running)

# Check service
sudo service mysql status
```

**Start/Stop MySQL if Needed**
```bash
# Start MySQL
sudo systemctl start mysql
sudo service mysql start

# Stop MySQL
sudo systemctl stop mysql

# Restart
sudo systemctl restart mysql

# Enable auto-start on boot
sudo systemctl enable mysql
```

#### Setup
**Found MySQL Credentials**
```bash
sudo cat /etc/mysql/debian.cnf
# DO NOT use debian-sys-maint for your application! This is an internal maintenance account that Debian/Ubuntu uses for system tasks.
```

**Setup Private User**
- Login with debian-sys-maint to set up root/user
```bash
mysql -u debian-sys-maint -p 
# Enter password obtained from sudo cat /etc/mysql/debian.cnf
```

- Set root password (if not set)
```sql
ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'secret';
FLUSH PRIVILEGES;

CREATE DATABASE IF NOT EXISTS kloe_bdt;

CREATE USER IF NOT EXISTS 'kloe_user'@'localhost' IDENTIFIED BY 'kloe_password';

GRANT ALL PRIVILEGES ON kloe_bdt.* TO 'kloe_user'@'localhost';

FLUSH PRIVILEGES;

SHOW DATABASES;
SELECT user, host FROM mysql.user;
SHOW GRANTS FOR 'kloe_user'@'localhost';

EXIT;
```

### Find the process using port 8000
```bash
sudo lsof -i :8000
```

# Kill all processes on port 8000
```bash
sudo fuser -k 8000/tcp
```