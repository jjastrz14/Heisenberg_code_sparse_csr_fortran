ifort -o helloworld_csr helloworld_csr.f90 -L/opt/FEAST/4.0/lib/x64 -lfeast -qmkl

ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90  main_spin_code.f90 -o heisenberg_main -L/opt/FEAST/4.0/lib/x64 -lfeast 
-qopenmp -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel -heap-arrays -integer-size 64 -check bounds


ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90  main_spin_code.f90 -o heisenberg_main -L/opt/FEAST/4.0/lib/x64 -lfeast -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel -heap-arrays -check bounds



to działa: ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90 main_spin_code.f90 -o heisenberg_main -L/opt/FEAST/4.0/lib/x64 -lfeast -qopenmp -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel 


source /opt/intel/oneapi/setvars.sh - for intel on macOS
