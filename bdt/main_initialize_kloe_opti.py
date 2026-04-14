import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)
import joblib
import uproot
import numpy as np
import pandas as pd
#from plots import plot_compr_hist, plot_var, plot_feature_pairs, plot_feature_target
from training_config import prepare_3photon_paris
from sklearn.model_selection import train_test_split
import random
import awkward as ak
import gc  # garbage collection

"""
python main_initialize_kloe_opti.py \
    --input /home/kloe/Desktop/KLOE_BDT/dataset/kloe_bdt.root \
    --chunk-size 50000 \
    --output-dir /home/kloe/Desktop/KLOE_BDT/dataset_bdt

OR

uv run main_initialize_kloe_opti.py \
    --input /home/kloe/Desktop/KLOE_BDT/dataset/kloe_bdt.root \
    --chunk-size 50000 \
    --output-dir /home/kloe/Desktop/KLOE_BDT/dataset_bdt
"""


def create_dataset(df, category): # For photon 4-momentum
    print(f'\n✅ Creating dataset ...')
    # Define photon 4-momentum
    br_nm = ['Br_E1', 'Br_px1', 'Br_py1', 'Br_pz1', 
             'Br_E2', 'Br_px2', 'Br_py2', 'Br_pz2', 
             'Br_E3', 'Br_px3', 'Br_py3', 'Br_pz3',
             'Br_m3pi', 'Br_lagvalue_min_7C', 'Br_deltaE',
             'Br_angle_pi0gam12', 'Br_ppIM', 'Br_betapi0',
             'Br_recon_indx', 'Br_bkg_indx']
    
    # Check which columns actually exist
    available_cols = [col for col in br_nm if col in df.columns]
    if len(available_cols) < len(br_nm):
        print(f"  Note: Using {len(available_cols)}/{len(br_nm)} available columns")
    
    # Selection cut, ensure physical region
    cut_region = (df['Br_lagvalue_min_7C'] < 100) if 'Br_lagvalue_min_7C' in df.columns else pd.Series(True, index=df.index)
    phys_region = ((df['Br_betapi0'] < 1) & (df['Br_betapi0'] > 0)) if 'Br_betapi0' in df.columns else pd.Series(True, index=df.index)

    #df = df[(df['Br_lagvalue_min_7C'] < 100) & (df['Br_betapi0'] < 1) & (df['Br_betapi0'] > 0)][br_nm]
    # Apply filters and create a proper copy
    df_filtered = df[cut_region & phys_region][available_cols].copy()
    df = df_filtered  # Now df is a clean copy

    # Create all_df, pos_df, neg_df for signal and background events
    if len(available_cols): # Check para length and br_nm length are consistent
        
        if category == 'signal':
            print(f"Creating all_df for {category} {df.columns}...")

            # ADDED: Check if classification columns exist
            if 'Br_recon_indx' in df.columns and 'Br_bkg_indx' in df.columns:
                pos_df = df[(df['Br_recon_indx'] == 2) & (df['Br_bkg_indx'] == 1)][available_cols]
                neg_df = df[~((df['Br_recon_indx'] == 2) & (df['Br_bkg_indx'] == 1))][available_cols]
            else:
                print(f"  Warning: Missing classification columns, treating all as background")
                pos_df = pd.DataFrame(columns=available_cols)
                neg_df = df[available_cols]

            # True positive
            nb_pos = [i for i in range(len(pos_df))]  
            if len(pos_df) > 0:
                pos_df = pos_df.copy()
                pos_df.insert(0, 'event', nb_pos)  # Add event column
                pos_df['is_signal'] = 1 # Add is_signal column
                pos_df['true_pi0_pair'] = [(0, 1)] * len(nb_pos)  # Add true_pi0_pair column
            
            # True negative
            nb_neg = [i for i in range(len(neg_df))]
            if len(neg_df) > 0:
                neg_df = neg_df.copy()
                neg_df.insert(0, 'event', nb_neg) # Add event column 
                neg_df['is_signal'] = 0 # Add is_signal column
                neg_df['true_pi0_pair'] = [(-1, -1)] * len(nb_neg) # Add true_pi0_pair column

            # Combine pos + neg dataset and shuffling
            dfs_to_concat = []
            if len(pos_df) > 0:
                dfs_to_concat.append(pos_df)
            if len(neg_df) > 0:
                dfs_to_concat.append(neg_df)
            
            if dfs_to_concat:
                all_df = pd.concat(dfs_to_concat, ignore_index=True)
                all_df = all_df.sample(frac=1).reset_index(drop=True)
            else:
                all_df = pd.DataFrame()
      
        elif category == 'background':
            print(f"Creating pho4mom_all_df for {category} {df.columns} ...")

            all_df = df.copy()
            nb_all_df = [i for i in range(len(all_df))]  
            all_df.insert(0, 'event', nb_all_df)  # Add event column
            all_df['is_signal'] = 0 # Add is_signal column
            all_df['true_pi0_pair'] = [(-1, -1)] * len(all_df) # Add true_pi0_pair column
        else: # combined or others
            raise ValueError("Only sig and bkg allow!")
    else: 
        print("No dataset other than signal or bkg or combined is expected!")
        raise ValueError("Array length mismatch - cannot proceed")

    # pi0 features for ML learning
    # MODIFIED: Only create pairs if we have data
    if len(all_df) > 0:
        pi0_all_df = prepare_3photon_paris(all_df)
    else:
        pi0_all_df = pd.DataFrame()
        
    return all_df, pi0_all_df

##
def data_splitting(all_df):
    print(f'\n✅ Splitting dataset ...')

    # First split: separate test set 20%
    all_df_trainval, all_df_test = train_test_split(
        all_df, test_size=0.2, random_state=42
    )

    # Second split: separte validation from training (20% of total)
    all_df_train, all_df_val = train_test_split(
        all_df_trainval, test_size=0.25, random_state=42
    )

    print(f"Training events:   {len(all_df_train)} ({len(all_df_train)/len(all_df)*100:.1f}%)")
    print(f"Validation events: {len(all_df_val)} ({len(all_df_val)/len(all_df)*100:.1f}%)")
    print(f"Test events:       {len(all_df_test)} ({len(all_df_test)/len(all_df)*100:.1f}%)")

    #print("CREATING PAIRS FROM EACH SPLIT")

    pair_train = prepare_3photon_paris(all_df_train)
    pair_val = prepare_3photon_paris(all_df_val)
    pair_test = prepare_3photon_paris(all_df_test)

    # Verify event column preservation
    for name, pair_df in [('Train', pair_train), ('Val', pair_val), ('Test', pair_test)]:
        if 'event' not in pair_df.columns:
            print(f"❌ CRITICAL: 'event' column missing in {name} pairs!")
            print(f"   Columns: {pair_df.columns.tolist()}")
        else:
            n_events = pair_df['event'].nunique()
            n_pairs = len(pair_df)
            print(f"✅ {name}: {n_pairs} pairs from {n_events} events ({n_pairs/n_events:.2f} pairs/event)")

    #print(f"Training pairs:   {len(pair_train)}")
    #print(f"Validation pairs: {len(pair_val)}")
    #print(f"Test pairs:       {len(pair_test)}")

    ## Step 3: Prepare features for training
    features = ['m_gg', 'opening_angle', 'cos_theta', 'E_asym', 'e_min_x_angle', 'E1', 'E2', 'E3', 'asym_x_angle', 'E_diff']
    
    # ADDED: Check which features actually exist
    available_features = [f for f in features if f in pair_train.columns]
    if len(available_features) < len(features):
        print(f"  Warning: Using {len(available_features)}/{len(features)} available features")
    
    X_train = pair_train[available_features] if len(pair_train) > 0 else pd.DataFrame()
    y_train = pair_train['is_pi0'] if 'is_pi0' in pair_train.columns else pd.Series()
    X_val = pair_val[available_features] if len(pair_val) > 0 else pd.DataFrame()
    y_val = pair_val['is_pi0'] if 'is_pi0' in pair_val.columns else pd.Series()
    X_test = pair_test[available_features] if len(pair_test) > 0 else pd.DataFrame()
    y_test = pair_test['is_pi0'] if 'is_pi0' in pair_test.columns else pd.Series()

    print(f"\n✅ Data ready for training:")
    print(f"  X_train shape: {X_train.shape}, y_train shape: {y_train.shape}")
    print(f"  X_val shape:   {X_val.shape}, y_val shape:   {y_val.shape}")
    print(f"  X_test shape:  {X_test.shape}, y_test shape:  {y_test.shape}")

    return all_df_train, all_df_val, all_df_test, X_train, y_train, X_val, y_val, X_test, y_test, pair_train, pair_val, pair_test

if __name__ == '__main__':
    #============================================================
    # LOAD INPUT ROOT FILES
    #============================================================
    # MODIFIED: Allow command line argument for input file
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default="../data/kloe_sample.root", 
                       help='Input ROOT file path')
    parser.add_argument('--chunk-size', type=int, default=50000, 
                       help='Number of entries to process at once')
    parser.add_argument('--max-events', type=int, default=None,
                       help='Maximum number of events to process (for testing)')
    parser.add_argument('--output-dir', type=str, default='./dataset_large',
                   help='Output directory for dataset files')
    args = parser.parse_args()
    
    f_nm = args.input

    # Create output directory
    data_dir = args.output_dir
    os.makedirs(data_dir, exist_ok=True)

    ## Loop over branches and create phys_map dynamically
    try:
        # Open the root file
        root_file = uproot.open(f_nm)
        branches = root_file.keys()
        print('All keys:', branches)
    
    except Exception as e:
        print('Error opening file: ', e)
        sys.exit(1)

    # Check first few braches and create a phys_map dynamically
    phys_map = {}

    for i, br_nm in enumerate(branches):
        # Remove ROOT cycle number (;1, ;2, etc.) for comparison
        base_br_nm = br_nm.split(';')[0]

        # Skip if already processed (in case of multiple cycles)
        if base_br_nm in phys_map:
            print(f"Skippin duplicate: {br_nm} (already have {base_br_nm})")
            continue

        #print(base_br_nm)
        if base_br_nm == "TISR3PI_SIG":
            br_title = rf"$e^{{+}}e^{{-}}\to\pi^{{+}}\pi^{{-}}\pi^{{0}}\gamma$"
            category = "signal"
        elif base_br_nm == "TOMEGAPI":
            br_title = rf"$\omega\pi^{0}$"
            category = "background"
        elif base_br_nm == "TKPM":
            br_title = rf"$e^{{+}}e^{{-}}\to\phi\to K\bar{{K}}$"
            category = "background"
        elif base_br_nm == "TKSL":
            br_title = rf"$e^{{+}}e^{{-}}\to\phi\to K_{{S}}K_{{L}}$"
            category = "background"
        elif base_br_nm == "TRHOPI":
            br_title = rf"$e^{{+}}e^{{-}}\to\phi\to \rho\pi$"
            category = "background"
        elif base_br_nm == "TBKGREST":
            br_title = rf"Others"
            category = "background"
        elif base_br_nm == "TEEG":
            br_title = rf"$e^{{+}}e^{{-}}\to\phi\to e^{{+}}e^{{-}}\gamma$"
            category = "background"
        elif base_br_nm == "TETAGAM":
            br_title = rf"$e^{{+}}e^{{-}}\to\phi\to\eta\gamma$"
            category = "signal"
        else:
            continue

        # Fill phys_map entries
        phys_map[base_br_nm] = {
            'br_title': br_title,
            'category': category,
            'original_name': br_nm # Store original for reference
        }
    
    # Check phys_map
    for data_type, info in phys_map.items():
        info_title = info['br_title']
        info_category = info['category']
        print(data_type, info_title, info_category)

    #============================================================
    # CREATE DATASET - MODIFIED to handle large files with chunking
    #============================================================
    df_list = []        # List for storing all dataset for combining
       

    df_train_list = []  # List for storing training splits
    df_val_list = []    # List for storing validation splits
    df_test_list = []   # List for storing test splits 
    pair_train_list = []
    pair_val_list = []
    pair_test_list = []
    X_train_list = []
    y_train_list = []
    X_val_list = []
    y_val_list = []
    X_test_list = []
    y_test_list = []

    ch_indx = 0
    for data_type, info in phys_map.items():
        try:
            data_nm = data_type.split(';')[0]

            ch_indx += 1
            print("="*25 + f"Channel {ch_indx}: {data_type}" + "="*25)
            
            # Get tree for this channel
            tree = root_file[data_type] 

            # Debug info for this channel
            print(f"Processing {data_type}...")
            total_entries = tree.num_entries
            if args.max_events:
                total_entries = min(total_entries, args.max_events)
            print(f"Number of entries: {total_entries}")

            # Read as awkward array to determine fields (only first 100)
            ak_array = tree.arrays(library="ak", entry_stop=100)

            # Generate fields_to_use for THIS channel
            exclude_fields = ['Br_pull_E1', 'Br_pull_x1', 'Br_pull_y1', 'Br_pull_z1', 'Br_pull_t1']
            fields_to_use = []

            for field in ak_array.fields:
                # Skip fields that match exclude_fields
                if any(pattern in field for pattern in exclude_fields):
                    continue

                # Only include 1D fields
                if ak_array[field].ndim == 1:
                    fields_to_use.append(field)
                else:
                    print(f"Excluding multi-dim field: {field} (ndim={ak_array[field].ndim})")

            # Check if we have the required branches
            required_br = ['Br_E1', 'Br_px1', 'Br_py1', 'Br_pz1', 
                        'Br_E2', 'Br_px2', 'Br_py2', 'Br_pz2', 
                        'Br_E3', 'Br_px3', 'Br_py3', 'Br_pz3',
                        'Br_m3pi', 'Br_lagvalue_min_7C', 'Br_deltaE',
                        'Br_angle_pi0gam12', 'Br_ppIM', 'Br_betapi0',
                        'Br_recon_indx', 'Br_bkg_indx']
            
            missing_br = [br for br in required_br if br not in fields_to_use]
            if missing_br:
                print(f"WARNING: Missing branches in {data_type}: {missing_br}")
                if info['category'] == 'signal' and len(missing_br) > 5:
                    print(f"Too many missing branches for signal channel {data_type}, skipping...")
                    continue
            
            print(f"Fields being used: {len(fields_to_use)} fields")

            # MODIFIED: Process in chunks for large files
            chunk_size = min(args.chunk_size, total_entries)
            n_chunks = (total_entries + chunk_size - 1) // chunk_size
            
            print(f"Processing {total_entries} entries in {n_chunks} chunks of {chunk_size}")
            
            all_chunks = []
            
            for chunk_idx in range(n_chunks):
                start = chunk_idx * chunk_size
                stop = min((chunk_idx + 1) * chunk_size, total_entries)
                
                print(f"\n  Processing chunk {chunk_idx+1}/{n_chunks} (entries {start}-{stop})")
                
                # Read chunk
                chunk_df = tree.arrays(fields_to_use, library="pd", 
                                      entry_start=start, entry_stop=stop)
                
                print(f"  Chunk shape: {chunk_df.shape}")
                
                # Process chunk
                chunk_all_df, chunk_pi0_df = create_dataset(chunk_df, info['category'])
                
                if chunk_all_df is not None and len(chunk_all_df) > 0:
                    # Add chunk identifier to event numbers to avoid duplicates
                    chunk_all_df['event'] = f"ch{chunk_idx}_" + chunk_all_df['event'].astype(str)
                    all_chunks.append(chunk_all_df)
                
                # Clear memory
                del chunk_df, chunk_all_df, chunk_pi0_df
                gc.collect()
            
            # Combine chunks
            if all_chunks:
                print(f"\nCombining {len(all_chunks)} chunks for {data_type}...")
                all_df = pd.concat(all_chunks, ignore_index=True)
                del all_chunks
                gc.collect()
                
                print(f"Total events after filtering: {len(all_df)}")
                
                # Check for anomalies
                if 'Br_betapi0' in all_df.columns:
                    betapi0_values = all_df['Br_betapi0']
                    print(f"Betapi0 stats: {betapi0_values.describe()}")
                
                # Data splitting
                if len(all_df) < 100:
                    print("WARNING: Very few events!")  
                    continue
                
                # MODIFIED: Create pi0 pairs after splitting (original behavior)
                all_df_train, all_df_val, all_df_test, X_train, y_train, X_val, y_val, X_test, y_test, pair_train, pair_val, pair_test = data_splitting(all_df)
                
                # Add to combined list
                df_list.append(all_df)
                df_train_list.append(all_df_train)
                df_val_list.append(all_df_val)
                df_test_list.append(all_df_test)

                pair_train_list.append(pair_train)
                pair_val_list.append(pair_val)
                pair_test_list.append(pair_test)
                X_train_list.append(X_train)
                y_train_list.append(y_train)
                X_val_list.append(X_val)
                y_val_list.append(y_val)
                X_test_list.append(X_test)
                y_test_list.append(y_test)

                # Save files with compression
                print(f"\nSaving files for {data_nm}...")
                
                joblib.dump(all_df, f'{data_dir}/all_df_{data_nm}.pkl', compress=3)
                
                if len(all_df_train) > 0:
                    joblib.dump(all_df_train, f'{data_dir}/all_df_train_{data_nm}.pkl', compress=3)
                if len(all_df_val) > 0:
                    joblib.dump(all_df_val, f'{data_dir}/all_df_val_{data_nm}.pkl', compress=3)
                if len(all_df_test) > 0:
                    joblib.dump(all_df_test, f'{data_dir}/all_df_test_{data_nm}.pkl', compress=3)
                
                if len(pair_train) > 0:
                    joblib.dump(pair_train, f'{data_dir}/pair_train_{data_nm}.pkl', compress=3)
                if len(pair_val) > 0:
                    joblib.dump(pair_val, f'{data_dir}/pair_val_{data_nm}.pkl', compress=3)
                if len(pair_test) > 0:
                    joblib.dump(pair_test, f'{data_dir}/pair_test_{data_nm}.pkl', compress=3)
                
                if len(X_train) > 0:
                    joblib.dump(X_train, f'{data_dir}/X_train_{data_nm}.pkl', compress=3)
                if len(X_val) > 0:
                    joblib.dump(X_val, f'{data_dir}/X_val_{data_nm}.pkl', compress=3)
                if len(X_test) > 0:
                    joblib.dump(X_test, f'{data_dir}/X_test_{data_nm}.pkl', compress=3)
                
                if len(y_train) > 0:
                    joblib.dump(y_train, f'{data_dir}/y_train_{data_nm}.pkl', compress=3)
                if len(y_val) > 0:
                    joblib.dump(y_val, f'{data_dir}/y_val_{data_nm}.pkl', compress=3)
                if len(y_test) > 0:
                    joblib.dump(y_test, f'{data_dir}/y_test_{data_nm}.pkl', compress=3)
                
                
            else:
                print(f"No valid events for {data_type}")
  
        except Exception as e:
            print(f"ERROR processing {data_type}: {str(e)}")
            import traceback
            traceback.print_exc()
            print(f"Skipping {data_type} and continuing with next channel...")
            continue

    ## Combining dataset
    if df_list:
        print(f"\nCombining {len(df_list)} channels ...")
        
        # Add channel prefix to event IDs to avoid collisions
        for i, (df, (data_type, info)) in enumerate(zip(df_list, list(phys_map.items())[:len(df_list)])):
            channel_name = data_type.split(';')[0]
            df['event'] = f"{channel_name}_" + df['event'].astype(str)
            print(f"    Updated event IDs for {channel_name}")

        # Combine full datasets
        df_comb = pd.concat(df_list, ignore_index=True)
        print(f"Raw combined shape: {df_comb.shape}")

        # Shuffle
        df_comb = df_comb.sample(frac=1, random_state=42).reset_index(drop=True)
        print(f"Shuffled combined shape: {df_comb.shape}")

        # Initialize combined split variables
        #all_df_train_comb = pd.DataFrame()
        #all_df_val_comb = pd.DataFrame()
        #all_df_test_comb = pd.DataFrame()
      

        # Combine training splits
        if df_train_list:
            print(f"\nCombining training splits from {len(df_train_list)} channels...")
            all_df_train_comb = pd.concat(df_train_list, ignore_index=True)
            all_df_train_comb = all_df_train_comb.sample(frac=1, random_state=42).reset_index(drop=True)
            print(f"Combined training events: {len(all_df_train_comb)}")

            # Combine training pairs if they exist
            if pair_train_list:
                pair_train_comb = pd.concat(pair_train_list, ignore_index=True)
                print(f"Combined training pairs: {len(pair_train_comb)}")
            
            # Combine training features and labels
            if X_train_list and y_train_list:
                X_train_comb = pd.concat(X_train_list, ignore_index=True)
                y_train_comb = pd.concat(y_train_list, ignore_index=True)
                print(f"X_train shape: {X_train_comb.shape}, y_train shape: {y_train_comb.shape}")
            
            # Save training split
            joblib.dump(all_df_train_comb, f'{data_dir}/all_df_train_TCOMB.pkl', compress=3)
            if pair_train_list:
                joblib.dump(pair_train_comb, f'{data_dir}/pair_train_TCOMB.pkl', compress=3)
            if X_train_list:
                joblib.dump(X_train_comb, f'{data_dir}/X_train_TCOMB.pkl', compress=3)
                joblib.dump(y_train_comb, f'{data_dir}/y_train_TCOMB.pkl', compress=3)

            
            # summary
            total_train = sum(len(df) for df in df_train_list)
            print(f"Training events:   {total_train}")
            # Print training info here where variables are defined
            if X_train_list and y_train_list:
                print(f"  Training pairs: {sum(len(pairs) for pairs in pair_train_list) if pair_train_list else 0}")
                print(f"  X_train components: {len(X_train_list)}")
            
        # Combine validation splits
        if df_val_list:
            print(f"\nCombining validation splits from {len(df_val_list)} channels...")
            all_df_val_comb = pd.concat(df_val_list, ignore_index=True)
            all_df_val_comb = all_df_val_comb.sample(frac=1, random_state=42).reset_index(drop=True)
            print(f"Combined validation events: {len(all_df_val_comb)}")

            if pair_val_list:
                pair_val_comb = pd.concat(pair_val_list, ignore_index=True)
                print(f"Combined validation pairs: {len(pair_val_comb)}")
        
            if X_val_list and y_val_list:
                X_val_comb = pd.concat(X_val_list, ignore_index=True)
                y_val_comb = pd.concat(y_val_list, ignore_index=True)
                print(f"X_val shape: {X_val_comb.shape}, y_val shape: {y_val_comb.shape}")

            joblib.dump(all_df_val_comb, f'{data_dir}/all_df_val_TCOMB.pkl', compress=3)

            if pair_val_list:
                joblib.dump(pair_val_comb, f'{data_dir}/pair_val_TCOMB.pkl', compress=3)

            if X_val_list:
                joblib.dump(X_val_comb, f'{data_dir}/X_val_TCOMB.pkl', compress=3)
                joblib.dump(y_val_comb, f'{data_dir}/y_val_TCOMB.pkl', compress=3)
            
            # summmary
            total_val = sum(len(df) for df in df_val_list)
            print(f"Validation events: {total_val}")
            if X_val_list and y_val_list:
                print(f"  Validation pairs: {sum(len(pairs) for pairs in pair_val_list) if pair_val_list else 0}")
                print(f"  X_val components: {len(X_val_list)}")
        
        # Combine test splits
        if df_test_list:
            print(f"\nCombining test splits from {len(df_test_list)} channels...")
            all_df_test_comb = pd.concat(df_test_list, ignore_index=True)
            all_df_test_comb = all_df_test_comb.sample(frac=1, random_state=42).reset_index(drop=True)
            print(f"Combined test events: {len(all_df_test_comb)}")

            if pair_test_list:
                pair_test_comb = pd.concat(pair_test_list, ignore_index=True)
                print(f"Combined test pairs: {len(pair_test_comb)}")
        
            if X_test_list and y_test_list:
                X_test_comb = pd.concat(X_test_list, ignore_index=True)
                y_test_comb = pd.concat(y_test_list, ignore_index=True)
                print(f"X_test shape: {X_test_comb.shape}, y_test shape: {y_test_comb.shape}")

            joblib.dump(all_df_test_comb, f'{data_dir}/all_df_test_TCOMB.pkl', compress=3)

            if pair_test_list:
                joblib.dump(pair_test_comb, f'{data_dir}/pair_test_TCOMB.pkl', compress=3)
        
            if X_test_list:
                joblib.dump(X_test_comb, f'{data_dir}/X_test_TCOMB.pkl', compress=3)
                joblib.dump(y_test_comb, f'{data_dir}/y_test_TCOMB.pkl', compress=3)

            total_test = sum(len(df) for df in df_test_list)
            print(f"Test events:       {total_test}")
            if X_test_list and y_test_list:
                print(f"  Test pairs: {sum(len(pairs) for pairs in pair_test_list) if pair_test_list else 0}")
                print(f"  X_test components: {len(X_test_list)}")

        # Create datasets
        all_df_comb = df_comb
        print("\nCreating pairs for combined dataset...")
        pi0_all_df_comb = prepare_3photon_paris(all_df_comb)

        # Split
        #all_df_train_comb, all_df_val_comb, all_df_test_comb, X_train_comb, y_train_comb, X_val_comb, y_val_comb, X_test_comb, y_test_comb, pair_train_comb, pair_val_comb, pair_test_comb = data_splitting(all_df_comb)

        print(f"\n✅ Data ready for training:")
    
    
    

        # Save combined files with compression
        print("\nSaving combined files...")
        
        joblib.dump(all_df_comb, f'{data_dir}/all_df_TCOMB.pkl', compress=3)
        joblib.dump(pi0_all_df_comb, f'{data_dir}/pi0_all_df_TCOMB.pkl', compress=3)

        #if len(all_df_train_comb) > 0:
        #    joblib.dump(all_df_train_comb, f'{data_dir}/all_df_train_TCOMB.pkl', compress=3)
        #if len(all_df_val_comb) > 0:
        #    joblib.dump(all_df_val_comb, f'{data_dir}/all_df_val_TCOMB.pkl', compress=3)
        #if len(all_df_test_comb) > 0:
            #joblib.dump(all_df_test_comb, f'{data_dir}/all_df_test_TCOMB.pkl', compress=3)

        #if len(pair_train_comb) > 0:
        #    joblib.dump(pair_train_comb, f'{data_dir}/pair_train_TCOMB.pkl', compress=3)
        #if len(pair_val_comb) > 0:
            #joblib.dump(pair_val_comb, f'{data_dir}/pair_val_TCOMB.pkl', compress=3)
        #if len(pair_test_comb) > 0:
        #    joblib.dump(pair_test_comb, f'{data_dir}/pair_test_TCOMB.pkl', compress=3)

        #if len(X_train_comb) > 0:
        #    joblib.dump(X_train_comb, f'{data_dir}/X_train_TCOMB.pkl', compress=3)
        #if len(X_val_comb) > 0:
        #    joblib.dump(X_val_comb, f'{data_dir}/X_val_TCOMB.pkl', compress=3)
        #if len(X_test_comb) > 0:
        #    joblib.dump(X_test_comb, f'{data_dir}/X_test_TCOMB.pkl', compress=3)

        #if len(y_train_comb) > 0:
        #    joblib.dump(y_train_comb, f'{data_dir}/y_train_TCOMB.pkl', compress=3)
        #if len(y_val_comb) > 0:
        #    joblib.dump(y_val_comb, f'{data_dir}/y_val_TCOMB.pkl', compress=3)
        #if len(y_test_comb) > 0:
        #   joblib.dump(y_test_comb, f'{data_dir}/y_test_TCOMB.pkl', compress=3)
        
        print(f"Combined data contains: {[k for k in phys_map.keys()]}")
        
        # Add TCOMB to phys_map
        phys_map['TCOMB'] = {
            'br_title': "MC combined",
            'category': "combined"
        }
    else:
        print(f"df_list is empty, no channels can be combined!")
        raise ValueError("No data to combine!")

    # Save phys_map
    joblib.dump(phys_map, f'{data_dir}/phys_map.pkl', compress=3)
    print(f"\n✅ All processing complete! Files saved to {data_dir}")
