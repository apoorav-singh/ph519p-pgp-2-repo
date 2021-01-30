# $ bin/bash
 #rm output.txt
 #gfortran -o output.o harmonic_wavefunction.f95
 gfortran -o output.o harmonic_numerov.f90
 ./output.o #>> output.txt
 python3 data_plot.py 