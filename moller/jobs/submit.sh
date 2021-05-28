#!/bin/bash


rm command1

touch command1

# rm *.out
python submit.py 1 11 >> command1
sh command1
sleep 10
