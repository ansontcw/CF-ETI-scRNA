#!/bin/bash

if [ $# -ne 2 ];then
        echo "Usage: sh $0 <capture> <workDIR>"
        exit 1
fi

capture=$1
workDIR=$2

module load cellranger/7.1.0
cd $workDIR
cellranger multi --id=${capture} --csv=/group/canc2/anson/working/CF-ETI-scRNA/data/sample_sheets/${capture}.config.csv --localcores=24 --localmem=64
