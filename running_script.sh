#$ bin/bash
 rm output.txt
 gfortran -o output.o harmonic_wavefunction.f90 
 ./output.o >> output.txt
 python3 data_plot.py 