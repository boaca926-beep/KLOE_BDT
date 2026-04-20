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

    

    
    
    

    