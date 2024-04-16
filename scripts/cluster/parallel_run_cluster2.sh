#!/usr/bin/env bash
# -*- encoding: utf-8 -*-

### Get the current date and start making the results directory
CURRDATE=$(date +"%FT%H%M")
RESDATA="/beegfs/nvardala/results_fig4_kN2/$1/data"
RESLOGS="/beegfs/nvardala/results_fig4_kN2/$1/logs/$1_$CURRDATE"

### Check and make directory for saving ALL the results
while : ; do
    if [ -d $RESLOGS ]; then
        echo "Directory $RESLOGS exists. Will add the results in here..."
        # sleep 5
        # continue
    else
        echo "Directory $RESLOGS does not exist. Creating..."
        mkdir $RESLOGS
        break
    fi
done

### Get all config directories
# CONF_DIRS=$(find ./configs/ -maxdepth 1 -mindepth 1 -type d)
# ISTIM=$1
# CONF_DIRS="configs/${ISTIM}_nA"
CONF_DIRS="configs/fig4_kN2/$1"

### Go through the config directories and do the following:
for DIR in $CONF_DIRS;
    do
        ### Store the basename for the parameters i.e. 10nA
        PNAME=$(basename $DIR)
        echo "In $PNAME"

        ### 1. Make a subdirectory to store the used configuration files for this parameter set i.e. results_{date}/10nA/...
        CURRRES="$RESLOGS/$PNAME"
        if ! [ -d $CURRRES ]; then
            echo "Directory $CURRRES does not exist. Creating..."
            mkdir $CURRRES
        fi

        BACKCONF="$CURRRES/configs"
        if ! [ -d $BACKCONF ]; then
            echo "Directory $BACKCONF does not exist. Creating..."
            mkdir $BACKCONF
        fi


        ### 2. Make the current config set output logs
        echo "Creating txt output for current config set..."
        touch "$CURRRES/total.txt"


        ### 3. Queue the batch script
        ### Set the variables used as arguments
        CNT=0
        RESDIR=$CURRRES
        for FN_CONF in $DIR/*.json
            do
                echo $FN_CONF

                ### Queue the simulation
                let DELAY="CNT*60"
            echo "Delayed: $DELAY seconds"
            DELAY_TOT=$(($DELAY+$2))
            echo "Delay TOTAL: $DELAY_TOT seconds"
                sbatch --job-name="SIM_${FN_CONF}" --begin=now+$DELAY_TOT --time=120 --cpus-per-task=1 --mem-per-cpu=24G --ntasks=1 --export=FCONF=$FN_CONF,RESDIR=$RESDIR,CNT=$CNT,OUTDIR=$RESDATA run_sim.sh

                let "CNT+=1"

                sleep 1

                ### Copy the config file used to results dir
                cp $FN_CONF $BACKCONF
            done
    done
