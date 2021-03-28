#!/bin/bash

Volumes="15. 25. 30. 60. 100. 200. 300."

for Vi in $Volumes; do
	wdir="w1k_int10k_140K_"$Vi
	mkdir $wdir
	cd $wdir

	../pypresso ../pvplotunit.py $Vi
	
	cd ../
done

exit 
