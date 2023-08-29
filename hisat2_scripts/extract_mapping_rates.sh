#!/bin/bash

echo -e "SRA \talignment_rate" > mapping_rates.tsv

for file in *.txt;
do
name=$(basename ${file});
last_line=$(tail -1 "$file" | tr -d '\n')
echo -e "${name%.*}\t${last_line}" >> mapping_rates.tsv
done

