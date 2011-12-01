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
    
    echo "MSA visualization.."
    flock $PROGRAMS_PATH/process_lock $PROGRAMS_PATH/msavisproject.py $psipredchoice $file_absolutepath $jobID $linewidth $PROGRAMS_PATH $format
    
    if [ ! $mail = "" ]
    then
      echo "sending email.."
      $PROGRAMS_PATH/mail.py $mail $PAGE_ADDRESS $jobID $datetime
    fi

    echo 'jobID' $jobID 'completed'
}

runvis
