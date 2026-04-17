#!/bin/bash
#INPATH=/media/bo/Backup/DATA_PEAK/ROOTINPUT #Backup
#INPATH=/media/bo/8E97-E8DD/DATA_PEAK/ROOTINPUT #8E97-E8DD
INPATH=/media/kloe/Analysis_Disk/ROOTINPUT #laptop kloe

OUTPATH=../path_norm

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

path_UFO1=$INPATH/UFO

echo "Creating input root list from" $INPATH
echo "outpath:" $OUTPATH

printf '%s\n' "$path_KSL1"/* "$path_KSL2"/* "$path_KSL3"/* "$path_KSL4"/* "$path_KSL5"/* "$path_KSL6"/* "$path_KSL7"/* "$path_KSL8"/* "$path_KSL9"/* > $OUTPATH/ksl_path

printf '%s\n' "$path_EXP1"/* "$path_EXP2"/* "$path_EXP3"/* "$path_EXP4"/* "$path_EXP5"/* "$path_EXP6"/* "$path_EXP7"/* "$path_EXP8"/* "$path_EXP9"/* > $OUTPATH/exp_path

printf '%s\n' "$path_EEG1"/* "$path_EEG2"/* "$path_EEG3"/* "$path_EEG4"/* "$path_EEG5"/* "$path_EEG6"/* "$path_EEG7"/* "$path_EEG8"/* "$path_EEG9"/* > $OUTPATH/eeg_path

printf '%s\n' "$path_SIG1"/* "$path_SIG2"/* "$path_SIG3"/* "$path_SIG4"/* "$path_SIG5"/* "$path_SIG6"/* "$path_SIG7"/* "$path_SIG8"/* > $OUTPATH/sig_path

printf '%s\n' "$path_UFO1"/* > $OUTPATH/ufo_path

DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")

for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    data_type=${DATA_TYPE[i]}
    
    total_lines=$(wc -l < $OUTPATH/${data_type}_path)  
    part_lines=$(( (total_lines + 2) / 3 ))  # Round up division
    
    echo ${data_type} $total_lines $part_lines
    
    # Part 1: Lines 1 to part_lines
    head -n $part_lines $OUTPATH/${data_type}_path > $OUTPATH/${data_type}_path1
    
    # Part 2: Lines (part_lines+1) to (2*part_lines)
    tail -n +$((part_lines + 1)) $OUTPATH/${data_type}_path | head -n $part_lines > $OUTPATH/${data_type}_path2
    
    # Part 3: Remaining lines
    tail -n +$((2 * part_lines + 1)) $OUTPATH/${data_type}_path > $OUTPATH/${data_type}_path3

done


