#!/bin/bash
for col_rows in 60 80 100
do
    for background in 30000 35000 40000 45000
    do
       python topology_2d.py $col_rows $col_rows $background "0.1" "0.1" "0.05" "0.075" "-5.0" "5.0" "-3000."
    done
done
