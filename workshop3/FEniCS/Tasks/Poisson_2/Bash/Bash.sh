#!/bin/bash
# Change filename variable and <problem>.py
declare -a nx=(5 10 20)

mv data.dat 0ld_data.dat
for i in "${nx[@]}"
do
echo "nx is: ${i}"

	filename="Possion_nx_${i}"
    	echo "Solve for $filename"
    	python ../code3.py ${i} > ./$filename.txt

    	awk '$1=="nxLOG" {print $3}' ./$filename.txt > ./results1_$filename.dat
    	awk '$1=="L2_errorLOG" {print $3}' ./$filename.txt > ./results2_$filename.dat

   	paste ./results1_$filename.dat ./results2_$filename.dat  >> ./data.dat
    	rm ./results1_$filename.dat ./results2_$filename.dat

	rm ./$filename.txt
done

#Plotting
python -c "import matplotlib.pyplot as plt;\
import numpy as np;data = np.loadtxt('data.dat');\
plt.plot(data[:,0], data[:,1], '-o');\
plt.xlabel('log(1/h-size)');plt.ylabel('log(L2 error)');\
plt.savefig('convergence.png');\
slope = (data[len(data)-1,1])-(data[0,1])-(data[len(data)-1,0]-data[0,0]);\
print slope " 

	
