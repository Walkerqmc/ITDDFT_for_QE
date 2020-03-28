# Imaginary-Time-Dependent-DFT-for-Quantum-ESPRESSO
A modification of the DFT program Quantum ESPRESSO that allows the user to propagate the wave function in imaginary time.  
For Ubuntu.

This modification gives the option of changing SCF iterations to imaginary time steps, where the wave function 
is propagated in imaginary time instead of running c_bands, which solves the energy eigenvalue equation. 

To setup:
  -first install Quantum ESPRESSO v6.3 https://github.com/QEF/q-e/releases/tag/qe-6.3 or use an existing version
  -download the modified files in this repository and place input_parameters.f90 in the folder MODULES 
   and the rest in the folder PW/src
  -open a terminal in the main folder and type './configure'  then 'make all'
  

Inputs for this modification are part of the "system" namelist on the input file, they are:

run_ITDDFT
   A logical value that will run ITDDFT if .true. and will run standard SCF if .false.
   
extensive_ITDDFT
   A logical value that if it and run_ITDDFT are .true. will super impose electrons among k-points.
   
g_dependant_dtau
   A logical value that if it and run_ITDDFT are .true. will propagate g-vectors that have less kinetic energy 
   with a greater time step.
   
dtau
   The time step length is equal to dtau multiplied by the inverse of the maximum g-vector kinetic energy.
   Except when g_dependant_dtau=.true., where The time step length is less than or equal to dtau multiplied
   by the inverse of the g-vector kinetic energy.
   
switch_tr
   A Threshold similar to conv_thr that when met will switch ITDDFT to SCF.
   
switch_iter
   Will switch ITDDFT to SCF after this many time steps.
   
freeze_band
   For extensive_ITDDFT=.true., it is the number of bands with full occupation number, starting from the lowest energy. 
