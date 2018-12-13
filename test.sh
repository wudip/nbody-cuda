#!/bin/bash

echo "input time"
for i in `seq 100000 100000 1000000`;
do
  ./generator $i > temp_input.csv
  TIME=`./nbody < temp_input.csv | tail -n 1 | cut -d" " -f 2`
  echo "$i $TIME"
done
rm temp_input.csv
