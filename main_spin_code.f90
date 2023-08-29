include "mkl_spblas.f90" !necessary for intel mkl blas for sparse matrix multiplicaton

program spin_code
    use spin_systems
    use tests_module
    use math_functions
    implicit none

    integer :: N_spin, N_spin_max, no_of_nonzero, size_of_sub_A, size_of_sub_B, i
    double precision :: J_spin 
    character(len = 12) :: N_spin_char, J_spin_char
    character(len=53) :: dec2bin, spin_basis
    integer, allocatable :: hash(:), indices_Sz_basis_sorted(:), target_sz(:), basis_vector(:,:), basis_rho_target(:,:), new_basis(:,:)
    double precision, allocatable :: Sz_basis(:)

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

    !call mmm_csr_test()
    call omp_mkl_small_test()
    call sparse_dfeast_test()
    call test_permutation_H_for_4_sites()

    write(*,*) '------- START Heisenberg Program -------'
    N_spin_max = 2**N_spin 
    allocate(hash(N_spin_max), target_sz(N_spin+1))
    allocate(indices_Sz_basis_sorted(N_spin_max), Sz_basis(N_spin_max), basis_vector(N_spin_max, N_spin))

    do i = 0, N_spin
        target_sz(i+1) = int(N_spin/2) - i
    end do 

    print *, "All combination of S_z spin", target_sz
    

    ! call H_create_basis_sz_sorted(N_spin, indices_Sz_basis_sorted)
    
    !creating a basis for N_spin system and whole basis vector -> 0000, 0001
    call H_create_basis_sz(N_spin, Sz_basis, basis_vector)
    
    ! here put loop for whole target_sz
    call Hash_basis_with_target(N_spin, target_sz(2), Sz_basis, basis_vector, hash, basis_rho_target)

    ! calulacting basis for subsystems A and B 

    call H_XXX_block_diag_with_target_dense(N_spin, J_spin, hash)
    call H_XXX_block_feast_vec_fill(N_spin, J_spin, hash, no_of_nonzero)
    call H_XXX_block_feast_vec_diag(N_spin, no_of_nonzero, hash)

    !spin_basis = dec2bin
    size_of_sub_A = int(N_spin/2)
    size_of_sub_B = int(N_spin - size_of_sub_A)
    print *, "Size of subsystem A", size_of_sub_A
    print *, "Size of subsystem B", size_of_sub_B

    call Basis_for_rho_reduced(N_spin, basis_rho_target, size_of_sub_A, size_of_sub_B, new_basis)

    deallocate(Sz_basis, basis_vector, new_basis)

    !call Basis_for_rho_reduced(spin_basis, size_of_sub_A, size_of_sub_B, new_basis)

    !steps: generuj macierze dla kazdego spinu i porównaj z pythonem - dziala 

    ! call CSR_matrix_multiplication_for_3_matrices(N_spin, J_spin, indices_Sz_basis_sorted)

    ! diagonalization using LAPACK from intel 
    !call H_XXX_diag(N_spin, J_spin)

    !call test_permutation_H_for_4_sites()

    ! diagonalization via FEAST algorithm (saving CSR matrices to file right now)
    !call sparse_dfeast_test()
    !call H_XXX_feast_vec_fill(N_spin, J_spin, no_of_nonzero)
    !call H_XXX_feast_vec_diag(N_spin, no_of_nonzero)

end program spin_code