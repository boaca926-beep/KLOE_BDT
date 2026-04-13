#!/bin/bash
#INPATH=/media/bo/Analysis_Disk/ROOTINPUT #Analysis_Disk
#INPATH=/media/bo/Backup/DATA_PEAK/ROOTINPUT #Backup
INPATH=/media/bo/8E97-E8DD/DATA_PEAK/ROOTINPUT #8E97-E8DD 
OUTPATH=../path_chain

# create the main folder
if [[ -d "$OUTPATH" ]]; then
    #echo "Remove ${OUTPATH} ..."
    rm -rf $OUTPATH
fi    

mkdir ${OUTPATH} # result folder

path_SIG1=$INPATH/SIG/isr3pi1
path_SIG2=$INPATH/SIG/isr3pi2
path_SIG3=$INPATH/SIG/isr3pi3
path_SIG4=$INPATH/SIG/isr3pi4
path_SIG5=$INPATH/SIG/isr3pi5
path_SIG6=$INPATH/SIG/isr3pi6
path_SIG7=$INPATH/SIG/isr3pi7
path_SIG8=$INPATH/SIG/isr3pi8

path_KSL1=$INPATH/KSL/mcksl1
path_KSL2=$INPATH/KSL/mcksl2
path_KSL3=$INPATH/KSL/mcksl3
path_KSL4=$INPATH/KSL/mcksl4
path_KSL5=$INPATH/KSL/mcksl5
path_KSL6=$INPATH/KSL/mcksl6
path_KSL7=$INPATH/KSL/mcksl7
path_KSL8=$INPATH/KSL/mcksl8
path_KSL9=$INPATH/KSL/mcksl9

path_EXP1=$INPATH/EXP/data1
path_EXP2=$INPATH/EXP/data2
path_EXP3=$INPATH/EXP/data3
path_EXP4=$INPATH/EXP/data4
path_EXP5=$INPATH/EXP/data5
path_EXP6=$INPATH/EXP/data6
path_EXP7=$INPATH/EXP/data7
path_EXP8=$INPATH/EXP/data8
path_EXP9=$INPATH/EXP/data9

path_EEG1=$INPATH/EEG/mceeg1
path_EEG2=$INPATH/EEG/mceeg2
path_EEG3=$INPATH/EEG/mceeg3
path_EEG4=$INPATH/EEG/mceeg4
path_EEG5=$INPATH/EEG/mceeg5
path_EEG6=$INPATH/EEG/mceeg6
path_EEG7=$INPATH/EEG/mceeg7
path_EEG8=$INPATH/EEG/mceeg8
path_EEG9=$INPATH/EEG/mceeg9

path_RHOPI1=$INPATH/RHOPI/MCRHOPI1
path_RHOPI2=$INPATH/RHOPI/MCRHOPI2
path_RHOPI3=$INPATH/RHOPI/MCRHOPI3
path_RHOPI4=$INPATH/RHOPI/MCRHOPI4
path_RHOPI5=$INPATH/RHOPI/MCRHOPI5
path_RHOPI6=$INPATH/RHOPI/MCRHOPI6
path_RHOPI7=$INPATH/RHOPI/MCRHOPI7
path_RHOPI8=$INPATH/RHOPI/MCRHOPI8
path_RHOPI9=$INPATH/RHOPI/MCRHOPI9

path_UFO1=$INPATH/UFO

echo "Creating input root list from" $INPATH
echo "outpath:" $OUTPATH

DATA_TYPE=("SIG" "KSL" "EXP" "EEG" "RHOPI")
MC_TYPE=("sig" "ksl" "exp" "eeg" "rhopi")
FILE_TYPE=("isr3pi1" "mcksl1" "data1" "mceeg1" "MCRHOPI1")

for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    data_type=${DATA_TYPE[i]}
    mc_type=${MC_TYPE[i]}
    file_type=${FILE_TYPE[i]}

    printf '%s\n' "$INPATH/${data_type}/${file_type}"/* > $OUTPATH/${mc_type}_path

    total_lines=$(wc -l < $OUTPATH/${mc_type}_path)  
    part_lines=$(( (total_lines + 2) / 3 ))  # Round up division
    
    echo $total_lines $part_lines
    
    # Part 1: Lines 1 to part_lines
    head -n $part_lines $OUTPATH/${mc_type}_path > $OUTPATH/${mc_type}_path1
    
    # Part 2: Lines (part_lines+1) to (2*part_lines)
    tail -n +$((part_lines + 1)) $OUTPATH/${mc_type}_path | head -n $part_lines > $OUTPATH/${mc_type}_path2
    
    # Part 3: Remaining lines
    tail -n +$((2 * part_lines + 1)) $OUTPATH/${mc_type}_path > $OUTPATH/${mc_type}_path3

done

printf '%s\n' "$path_UFO1"/* > $OUTPATH/ufo_path

total_lines=$(wc -l < $OUTPATH/ufo_path)  
part_lines=$(( (total_lines + 2) / 3 ))  # Round up division

echo $total_lines $part_lines

# Part 1: Lines 1 to part_lines
head -n $part_lines $OUTPATH/ufo_path > $OUTPATH/ufo_path1

# Part 2: Lines (part_lines+1) to (2*part_lines)
tail -n +$((part_lines + 1)) $OUTPATH/ufo_path | head -n $part_lines > $OUTPATH/ufo_path2

# Part 3: Remaining lines
tail -n +$((2 * part_lines + 1)) $OUTPATH/ufo_path > $OUTPATH/ufo_path3
