#!/bin/bash


rm command100
rm command200
rm command300
rm command400
rm command500
rm command600
rm command900
rm command1000
rm command1100



touch command100
touch command200
touch command300
touch command400
touch command500
touch command600
touch command700
touch command800
touch command900
touch command1000
touch command1100


# rm *.out
python submit.py "" "" "" 1 101 >> command100
sh command100
sleep 300

# rm *.out
python submit.py "" "" "" 101 201 >> command200
sh command200
sleep 300

# rm *.out
python submit.py "" "" "" 201 301 >> command300
sh command300
sleep 300

# rm *.out
python submit.py "" "" "" 301 401 >> command400
sh command400
sleep 300

# rm *.out
python submit.py "" "" "" 401 501 >> command500
sh command500
sleep 300


# rm *.out
python submit.py "" "" "" 501 601 >> command600
sh command600
sleep 300

# rm *.out
python submit.py "" "" "" 601 701 >> command700
sh command700
sleep 300

# rm *.out
python submit.py "" "" "" 701 801 >> command800
sh command800
sleep 300

# rm *.out
python submit.py "" "" "" 801 901 >> command900
sh command900
sleep 300

# rm *.out
python submit.py "" "" "" 901 1001 >> command1000
sh command1000
sleep 300

rm *.out
python submit.py "" "" "" 1001 1101 >> command1100
sh command1100
sleep 300





