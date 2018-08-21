#!/bin/bash
for col_rows in 80 100
do
    for background in 30000
    do
        #Columns,Rows, Background Rate, Radius excitational, Radius inhibitory, Sigma excitational, Sigma inhibitory, Inhibitory Weight, Excitational Weight
       python topology_2d.py $col_rows $col_rows $background "0.1" "0.1" "0.05" "0.075" "-5.0" "5.0" "-3000."
       python topology_2d.py $col_rows $col_rows $background "0.2" "0.2" "0.1" "0.1" "-5.0" "5.0" "-3000."
       python topology_2d.py $col_rows $col_rows $background "0.1" "0.2" "0.05" "0.1" "-5.0" "5.0" "-3000."
       python topology_2d.py $col_rows $col_rows $background "0.2" "0.1" "0.1" "0.05" "-10.0" "10.0" "-3000."
       python topology_2d.py $col_rows $col_rows $background "0.1" "0.2" "0.5" "0.75" "-20.0" "10.0" "-3000."
    done
done
