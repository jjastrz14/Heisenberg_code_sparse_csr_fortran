include "mkl_spblas.f90" !necessary for intel mkl blas for sparse matrix multiplicaton

program spin_code
    use spin_systems
    use tests_module
    use math_functions
    implicit none

    integer :: N_spin, N_spin_max, no_of_nonzero
    double precision :: J_spin 
    character(len = 12) :: N_spin_char, J_spin_char
    integer, allocatable :: hash(:), indices_Sz_basis_sorted(:), target_sz(:)

    ! Version 20.07 - dodajemy tworzenie bazy i liczenie entropii
    ! tworzenie bazy i s_z działa
    ! macierz permutacji H działa 
    ! test dla obrotu H n = 4 działa 

    ! do zrobienia
    ! zapis eigenvectora do pliku - podział na podpliki 
    
    ! TO DO LATER
    ! Zacznij ciac odpowiednie wektory do S_z , pytanie czy to w fortranie czy nie? 
    ! chyba w fotranie + diagonalizacja za pomocą feasta blockowych macierzy 
    ! w pythonie robisz to tnac pełną macierz H
    ! pytanie do Pawła: jak brać wycinek macierzy w H blockowym do feasta (numer eigenvalues i okno energetyczne)
    ! net cdf - binarna baza danych
   
    If(command_argument_count().NE.2) Then
        write(*,*)'Error, Only N (integer) and J (double precision) is required, program stopped'
        stop 
    endif 

    call get_command_argument(1, N_spin_char)
    call get_command_argument(2, J_spin_char)

    read(N_spin_char, *) N_spin
    read(J_spin_char, *) J_spin

    !N_spin = 4
    !J_spin = 1 

    write(*,*) 'Calculation of Heisenberg chain for N =', N_spin , 'and J = ', J_spin 
    write(*,*) ' '

    call mmm_csr_test()
    call omp_mkl_small_test()
    call sparse_dfeast_test()
    call test_permutation_H_for_4_sites()

    write(*,*) '---- START Heisenberg Program ----'
    N_spin_max = 2**N_spin 
    allocate(hash(N_spin_max), indices_Sz_basis_sorted(N_spin_max), target_sz(5) )

    target_sz = [2, 1, 0, -1, -2]

    ! call H_create_basis_sz(N_spin, indices_Sz_basis_sorted)

    call H_create_basis_sz_with_target(N_spin, hash, target_sz(3))
    call H_XXX_diag_with_target_dense(N_spin, J_spin, hash)

    !steps: generuj macierze dla kazdego spinu i porównaj z pythonem 

    ! call CSR_matrix_multiplication_for_3_matrices(N_spin, J_spin, indices_Sz_basis_sorted)

    ! diagonalization using LAPACK from intel 
    !call H_XXX_diag(N_spin, J_spin)

    !call test_permutation_H_for_4_sites()

    ! diagonalization via FEAST algorithm (saving CSR matrices to file right now)
    !call sparse_dfeast_test()
    !call H_XXX_feast_vec_fill(N_spin, J_spin, no_of_nonzero)
    !call H_XXX_feast_vec_diag(N_spin, no_of_nonzero)

end program spin_code