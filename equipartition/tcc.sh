#!/usr/bin/bash

# make a table called time.dat, which is a list of time corresponding snapshots in nbody units

#mkdir sigma
#cd sigma;mkdir 1 2 3 4 5 6 7 8 9 10
#cd ../

for name in $(ls .. --hide=scripts)
do
	cd ../$name
	rm -f *.R *.Rout time.dat filename.dat plottime.dat

#find time
	ls c_????.dat > filename.dat 
	for time in $(ls c_????.dat)
	do
		head -n 1 $time | tr '=' '\n' | tr '.' '\n' | tr -d " " | tail -n 2 | head -n 1 >> time.dat
	done

done
