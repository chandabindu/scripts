#!/bin/bash


rm command100

touch command100

# rm *.out
python submit.py 1 11 >> command100
sh command100
sleep 10
