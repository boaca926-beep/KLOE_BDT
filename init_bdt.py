from mysql_db import MySQLKLOEDB
import mysql.connector # mysql-connector-python
from mysql.connector import Error
import pandas as pd

def load_model_metadata(db):
    """Load model metadata into database"""

    import joblib
    import os
    from datetime import datetime

    # Path to the ML model
    model_path = "models/pi0_classifier_model_TCOMB.pkl"
    print(f"\n\nℹ️  Loading the model '{model_path}'")

    if not os.path.exists(model_path):
        print(f"⚠️  Model not found at {model_path}")
        print("   Skipping metadata insertion")
        return

    # Error handling
    try:
        # Loac the model to extract info
        model = joblib.load(model_path)

        # Get model info
        n_features = model.n_features_in_ if hasattr(model, 'n_features_in_') else 0 # Number of features
        feature_names = model.feature_names_in_.tolist() if hasattr(model, 'feature_names_in_') else []
        print(f"{n_features} features:\n {feature_names}")
        
        # Check if model is XGBoost
        model_type = type(model).__name__

        # Get model performance metric if it exists
        auc_score = None
        accuracy = None
        threshold = 0,5 # Default threshold

        # Load metric file
        metrics_path = "models/metrics_TCOMB.json"
        if os.path.exists(metrics_path):
            import json
            with open(metrics_path, 'r') as f:
                metrics = json.load(f)
                auc_score = metrics.get('auc')
                accuracy = metrics.get('accuracy')

        #print(f"accuracy: {accuracy}; auc: {auc_score}")
        
        db.conn.commit()
        print(f"✅ Model metadata inserted:")
        print(f"    - Model: pi0_classifier v1.0")
        print(f"    - Type: {model_type}")
        if auc_score:
            print(f"    - AUC: {auc_score:.4f}")
        

    except Exception as e:
        print(f"❌ Error loading model metadata: {e}")
    

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
                        # Show table structure
                        cursor.execute(f"SELECT * FROM {db[0]}.{table_name}")
                        content = cursor.fetchall()
                        if content:
                            # Get column names
                            cursor.execute(f"DESCRIBE {db[0]}.{table_name}")
                            columns = cursor.fetchall()
                            col_names = [col[0] for col in columns]

                            # Print column headers
                            print(f"\t\t Columns: {', '.join(col_names)}")

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
        confirm = input("Type 'YES' to initialize 'kloe_bdt' database : ")

        if confirm == "YES":

            # Initialize database with explicit parameters
            print("\n\nℹ️  INITIALIZING KLOE_BDT DATABASE...")
            with MySQLKLOEDB(
                host='localhost',
                user='kloe_user',
                password='kloe_password'
            ) as db:
                # Create database and tables
                db.create_database()
                db.create_tables()

                # Insert model metadata
                load_model_metadata(db)

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
                    },
                    # ========== TWO UNLABELED DATA POINTS (NULL) ==========
                    {
                        'event': (12346, 10001, 0.45, None),  # NULL/Unlabeled
                        'photon_pair': {
                            'invariant_mass': 140.5,
                            'opening_angle': 0.65,
                            'energy_asymmetry': 0.25,
                            'energy_ratio': 0.75,
                            'energy_difference': 25.3,
                            'min_energy_angle': 55.1,
                            'asymmetry_angle': 0.22,
                            'bdt_prediction': 0.45
                        }
                    },
                    {
                        'event': (12346, 10002, 0.38, None),  # NULL/Unlabeled
                        'photon_pair': {
                            'invariant_mass': 142.1,
                            'opening_angle': 0.58,
                            'energy_asymmetry': 0.31,
                            'energy_ratio': 0.69,
                            'energy_difference': 30.8,
                            'min_energy_angle': 62.3,
                            'asymmetry_angle': 0.28,
                            'bdt_prediction': 0.38
                        }
                    }
                ]
                
                print(f"\n\nℹ️  INSERTING TEST EVENTS...")
                for data in events_with_pairs:
                    event_id = db.insert_event(*data['event'])
                    db.insert_photon_pair(event_id, data['photon_pair'])
                    print(f"  Inserted event {event_id} with BDT score {data['event'][2]}")
                
                # See all events
                all_events = pd.read_sql("SELECT * FROM events", db.conn)
                print("\n\t\t=== All Events ===")
                print(all_events.to_string())

                # See all photon pairs
                all_pairs = pd.read_sql("SELECT * FROM photon_pairs", db.conn)
                print("\n\t\t=== All Photon Pairs ===")
                print(all_pairs.to_string())

                # Join query for high score events
                bdt_score = 0.9
                print(f"\n\t\t=== SELECTED EVENTS WITH BDT_SCORE > {bdt_score} ===")
                results = db.query_signal_events(min_score=0.9)
                if not results.empty:
                    print(results.to_string())
                else:
                    print("No events found with BDT score > 0.9")

                # Test get_training_data
                print(f"\n\nℹ️  TESTING GET TRAINING DATA (EXCLUDE NULL SIGNAL TYPE) ...")
                training_data = db.get_training_data(limit=100)
                print(f"Retrieved {len(training_data)} training events passed the selection")
                if not training_data.empty:
                    print("\nTraining data sample:")
                    print(training_data.head())
                    print(f"\nLabel distribution:")
                    print(training_data['label'].value_counts())
                else:
                    print("No training data found (only unlabeled events?)")

        else:
            print("❌ Quit")