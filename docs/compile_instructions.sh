ifort -o helloworld_csr helloworld_csr.f90 -L/opt/FEAST/4.0/lib/x64 -lfeast -qmkl

ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90  main_spin_code.f90 -o heisenberg_main -L/opt/FEAST/4.0/lib/x64 -lfeast 
-qopenmp -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel -heap-arrays -integer-size 64 -check bounds


ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90  main_spin_code.f90 -o heisenberg_main -L/opt/FEAST/4.0/lib/x64 -lfeast -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel -heap-arrays -check bounds



#to działa: 
ifort math_func_module.f90 spin_systems_module.f90 tests_module.f90 main_spin_code.f90 -o heisenberg_main -L/opt/FEAST/4.0/lib/x64 -lfeast -qopenmp -qmkl=cluster -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel 


source /opt/intel/oneapi/setvars.sh - for intel on macOS

#godot:

ifx -o heisenberg -module build/ -L$FEASTROOT -lfeast -qmkl -qopenmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -xHost -O3 -static-intel -heap-arrays -integer-size 64 -check bounds src/main_spin_code_feast_time_window.f90

#godot shrinked version: 
ifx -o heisenberg src/main_spin_code_feast_time_window.f90  -module build/ -L$FEASTROOT -lfeast -qmkl -qopenmp

#godot multi modules version: 

ifx -o heisenberg -module build/ src/utils.f90 src/tests_module.f90 src/spin.f90 src/main.f90 -L$FEASTROOT -lfeast -qmkl -qopenmp

#godot to run on several processors: taskset -c 0-3 ./heisenberg 10 1.0


mpirun -np 2 -genv OMP_NUM_THREADS 8 ./heisenberg_pfeast 4 1.0