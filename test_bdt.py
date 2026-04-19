from mysql_db import MySQLKLOEDB
import mysql.connector # mysql-connector-python
from mysql.connector import Error

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

            check_database_exists()
        else:
            print("Operation cancelled")
    
    else:
        #print("Database 'kloe_bdt' doesn't exist!")
        confirm = input ("Type 'YES' to initialize 'kloe_bdt' database : ")

        if confirm == "YES":

            # initialize database
            with MySQLKLOEDB() as db:
                db.create_database()
                db.create_tables()

                # Insert data
                event_id = db.insert_event(12345, 67890, 0.95, True)
                db.insert_photon_pair(event_id,{
                    'invariant_mass': 135.2
                })

                check_database_exists()

        else:
            print("❌ Quit")

    

    
    
    

    