from mysql_db import MySQLKLOEDB
import mysql.connector
from mysql.connector import Error
import pandas as pd
import os  # ADD THIS

# Database connection parameters - support both local and Docker
DB_HOST = os.getenv('MYSQL_HOST', 'localhost')
DB_USER = os.getenv('MYSQL_USER', 'kloe_user')
DB_PASSWORD = os.getenv('MYSQL_PASSWORD', 'kloe_password')
DB_ROOT_PASSWORD = os.getenv('MYSQL_ROOT_PASSWORD', 'secret')


def load_model_metadata(db):
    """Load model metadata into database"""
    import joblib
    import os as os_module
    import json
    from datetime import datetime

    # Path to the ML model
    model_path = "models/pi0_classifier_model_TCOMB.pkl"
    print(f"\n\nℹ️  Loading the model '{model_path}'")

    if not os_module.path.exists(model_path):
        print(f"⚠️  Model not found at {model_path}")
        print("   Skipping metadata insertion")
        return

    try:
        # Load the model to extract info (FIXED: removed duplicate with block)
        model = joblib.load(model_path)

        # Get model info
        n_features = model.n_features_in_ if hasattr(model, 'n_features_in_') else 0
        feature_names = []
        if hasattr(model, 'feature_names_in_'):
            feature_names = model.feature_names_in_.tolist()
        training_time = getattr(model, 'training_time_minutes_in_', 0)
        model_type = type(model).__name__

        # Get model performance metrics
        auc_score = None
        accuracy = None
        gpu_enabled = False
        threshold = 0.5

        # Load metric file
        metrics_path = "models/metrics_TCOMB.json"
        if os_module.path.exists(metrics_path):
            with open(metrics_path, 'r') as f:
                metrics = json.load(f)
                auc_score = metrics.get('auc')
                accuracy = metrics.get('accuracy')
                gpu_enabled = metrics.get('gpu_enabled', False)
                metrics_feature_names = metrics.get('feature_names', [])
                if training_time == 0:
                    training_time = metrics.get('training_time_minutes', 0)

                # Use feature names from metrics if not in model
                if not feature_names and metrics_feature_names:
                    feature_names = metrics_feature_names

                print(f"    Loaded metrics from {metrics_path}")
                print(f"    All keys: {list(metrics.keys())}")
                print(f"    GPU from metrics: {gpu_enabled}")
                print(f"    Training time: {training_time:.2f} minutes")
        else:
            print(f"    ⚠️  Metrics file not found at {metrics_path}")

        # Check if metadata already exists
        cursor = db.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM model_metadata WHERE model_name = 'pi0_classifier'")
        count = cursor.fetchone()[0]

        if count > 0:
            print("ℹ️  Updating existing metadata ...")
            cursor.execute("""
                UPDATE model_metadata
                SET model_version = 'v1.1',
                    n_features = %s,
                    feature_names = %s,
                    auc_score = %s,
                    accuracy = %s,
                    gpu_enabled = %s,
                    training_time_minutes = %s,
                    updated_at = NOW()
                WHERE model_name = 'pi0_classifier'
            """, (n_features, ','.join(feature_names), auc_score, accuracy, 
                  gpu_enabled, training_time))
        else:
            print("ℹ️  Inserting new metadata ...")
            cursor.execute("""
                INSERT INTO model_metadata (
                    model_name, model_version, is_active, training_date,
                    n_features, feature_names, auc_score, accuracy,
                    gpu_enabled, training_time_minutes
                ) VALUES (%s, %s, %s, NOW(), %s, %s, %s, %s, %s, %s)
            """, ('pi0_classifier', 'v1.0', True, n_features,
                  ','.join(feature_names), auc_score, accuracy,
                  gpu_enabled, training_time))

        db.conn.commit()
        print(f"\n✅ Model metadata inserted:")
        print(f"    - Model: pi0_classifier v1.0")
        print(f"    - Type: {model_type}")
        print(f"    - GPU enabled: {gpu_enabled}")

        if feature_names:
            preview = ', '.join(feature_names[:5])
            if len(feature_names) > 5:
                preview += "..."
            print(f"    - Features: {preview} ({len(feature_names)} total)")
        if auc_score:
            print(f"    - AUC: {auc_score:.4f}")
        if accuracy:
            print(f"    - Accuracy: {accuracy:.4f}")

    except Exception as e:
        print(f"❌ Error loading model metadata: {e}")
        import traceback
        traceback.print_exc()


def check_database_exists():
    """Check if kloe_bdt database exists"""
    try:
        conn = mysql.connector.connect(
            host=DB_HOST,
            user=DB_USER,
            password=DB_PASSWORD
        )
        cursor = conn.cursor()

        cursor.execute("SHOW DATABASES")
        databases = cursor.fetchall()

        for db in databases:
            if db[0] == 'kloe_bdt':
                print(f"Database {db[0]} EXISTS")
                cursor.execute(f"SHOW TABLES FROM {db[0]}")
                tables = cursor.fetchall()

                if tables:
                    print(f"Tables in {db[0]}")
                    for table in tables:
                        table_name = table[0]
                        print(f"\t- {table_name}")
                        cursor.execute(f"DESCRIBE {db[0]}.{table_name}")
                        columns = cursor.fetchall()
                        col_names = [col[0] for col in columns]
                        print(f"\t\t Columns: {', '.join(col_names)}")

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
    """Drop database using root (more privileges)"""
    try:
        conn = mysql.connector.connect(
            host=DB_HOST,
            user='root',
            password=DB_ROOT_PASSWORD
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
        print("\n⚠️  All data in kloe_bdt will be deleted!")
        confirm = input("Type 'DROP' to delete: ")

        if confirm == 'DROP':
            print("Dropping 'kloe_bdt' database")
            drop_database_via_root()
        else:
            print("Operation cancelled")

        check_database_exists()

    else:
        confirm = input("Type 'YES' to initialize 'kloe_bdt' database: ")

        if confirm == "YES":
            print("\n\nℹ️  INITIALIZING KLOE_BDT DATABASE...")
            with MySQLKLOEDB(
                host=DB_HOST,
                user=DB_USER,
                password=DB_PASSWORD
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
                    # Unlabeled data points (signal_type = NULL)
                    {
                        'event': (12346, 10001, 0.45, None),
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
                        'event': (12346, 10002, 0.38, None),
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
                print(f"\n\t\t=== SELECTED EVENTS WITH BDT_SCORE > 0.9 ===")
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