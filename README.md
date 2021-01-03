#  Project Title: ph519p-pgp-2-repo
This repository cantains part of the code for simulation program that have been part of my academic course work (PH-519P).

# Description
This branch of code cantains the code for the Numerical solution of Quantum Harmonic Oscillator using Numerov's Method. The code file "harmonic_numerov.f90" is taken from the lecture notes of Paolo Giannozzi (University of Udine). The second file "harmonic_wavefunction.f95" is being written by me.

## Running this code 
### harmonic_numerov.f90
There is Bash script file which will run the code for the case of harmonic_numerov.f90 comment out the output to output.txt file and subsequently the compilation of the other file (harmonic_wavefunction.f95) as it has inbuilt command which will generate the ouput file.

### harmonic_wavefunction.f95
As the case above same bash script file will run the code for the case of harmonic_wavefunction.f95 un-comment the output to output.txt command if commented. Then comment second fortran file (harmonic_numerov.f90), we don't want that both files run together. Output will be generted in output.txt file.

### Plotting
There are two option one is made available in the script i.e. using python and second is using Matlab for which I have included the file in the branch.