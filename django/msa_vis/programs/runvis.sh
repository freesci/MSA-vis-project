#!/bin/bash

psipredchoice=$1
file_absolutepath=$2
jobID=$3
linewidth=$4
PROGRAMS_PATH=$5
mail=$6
PAGE_ADDRESS=$7
datetime=$8
format=$9


function runvis()
{
    echo "lock acquire"
    $PROGRAMS_PATH/simple_lock.py acquire $PROGRAMS_PATH
    #exec 200<process_lock
    
    echo "MSA visualization.."
    $PROGRAMS_PATH/msavisproject.py $psipredchoice $file_absolutepath $jobID $linewidth $PROGRAMS_PATH $format
    
    if [ ! $mail = "" ]
    then
      echo "sending email.."
      $PROGRAMS_PATH/mail.py $mail $PAGE_ADDRESS $jobID $datetime
    fi

    echo "lock release"
    $PROGRAMS_PATH/simple_lock.py release $PROGRAMS_PATH
    #flock -n 200 || exit 1

    echo 'jobID' $jobID 'completed'
}

runvis
