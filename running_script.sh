#$ bin/bash

rm output.txt
gfortran -o output.o main_code.f90
./output.o >> output.txt
python3 data_plot.py 