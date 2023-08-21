# Heisenberg spin models and entalgement entropy measurments 
# Fotran version made with Maciej Bieniek

## Research Project: Twistronics - research on new quantum simulators

## compilation by ifort: 
*ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90  main_spin_code.f90 -o heisenberg_main -qopenmp -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel -heap-arrays -integer-size 64*

Field: Magnetism and strongly correlated systems

This repository contains work done in the pending research project at Nicolaus Copernicus University (UMK) in Toruń. 

List of topics:  
* 1D chain Heisenberg spin generator
* Graphs of Heisenberg spins 
* Density matrix from Heisenberg spin chain 
* Reduced density matrix calculation from Heisenberg 1D spin chain 
* Calculation of entalgement entropy of the corresponding subsytems
* Generation of Heisenberg Hamiltonian for CSR3 format for spare matrices
* Diagonalization of spare matrices by different methods: Lanczos and FEAST
*

JJ 2022