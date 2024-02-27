include "mkl_spblas.f90" !necessary for intel mkl blas for sparse matrix multiplicaton

program spin_code
    use spin_systems
    use tests_module
    use math_functions
    implicit none

    integer :: N_spin, N_spin_max, no_of_nonzero, size_of_sub_A, size_of_sub_B, i, k, target_sz_spin, number_of_eigen, Sz_subspace_size
    double precision :: J_spin, entropy_value, eigen_value_check, e_up, e_down
    double precision :: start, finish, finish_one_spin, start_one_spin, Sz_choice
    character(len = 12) :: N_spin_char, J_spin_char
    character(len=53) :: dec2bin, spin_basis, file_name1
    integer, allocatable :: Sz_basis(:), hash_Sz(:), hash(:), indices_Sz_basis_sorted(:), basis_vector(:,:), basis_rho_target(:,:), new_basis(:,:)
    double precision, allocatable :: target_sz(:), eigen_values(:), eigen_vectors(:,:), rho_reduced(:,:)

    ! OPIS DO POPRAWY! 
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

    !19.12 - wprowadzamy openMP
    !od teraz wersja kodu Maćka 
    ! entropie liczymy w pythonie 


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

    ! call mmm_csr_test()
    call omp_mkl_small_test()
    ! call test_permutation_H_for_4_sites()
    call sparse_zfeast_test()
    call sparse_dfeast_test()

    call cpu_time(start)
    write(*,*) '------- START Heisenberg Program -------'
    N_spin_max = 2**N_spin 
    allocate(hash_Sz(N_spin_max), target_sz(N_spin+1))
    allocate(indices_Sz_basis_sorted(N_spin_max), Sz_basis(N_spin_max), basis_vector(N_spin_max, N_spin))

    !write(*,*) '------- CSR Heisenberg Program -------'
    !call H_create_basis_sz_sorted(N_spin, indices_Sz_basis_sorted)
    !call CSR_matrix_multiplication_for_3_matrices(N_spin, J_spin, indices_Sz_basis_sorted)
    !write(*,*) '------- CSR END Heisenberg Program -------'

    ! writing all combination of Sz in decreasing order, works only for integers
    ! if float S_z is expected change it!
    do i = 0, N_spin
        target_sz(i+1) = N_spin/2 - i
    end do 

    print *, "All combination of S_z spin", target_sz
    write(*,*) " "

    file_name1 = 'Eigenvalues_results_' // trim(adjustl(N_spin_char)) // '_feast.dat'
    open (unit=1, file= trim(file_name1), recl=512)
    write(1,*), "# Eigenvalue     Spin_z    Sum_of_lambdas"

    Sz_choice = 0.0d0 ! integer counted from Sz_max (in a sense that Sz_choice = 1 means eg for N_spin=4, Sz_max=2 Sz=2)
    call Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)
    ! above we confirmed Sz_choice works

    call H_XXX_subspace_fill_and_diag_fast(N_spin, J_spin, Sz_subspace_size, hash_Sz) !here only csr3

    !do k = 1, size(target_sz, 1) ! loop over all S_z basis
        !chosen target in this loop
        !call cpu_time(start_one_spin)
        !target_sz_spin = target_sz(k)

        !creating a hashed basis for N_spin system and whole basis vector -> 0000, 0001 with specified target_sz
        !call Sz_subspace_choice(N_spin, target_sz[i], hash_Sz, Sz_subspace_size)
        !print *, "Basis created"
        !write(*,*) " "

        !write(*,*) 'Summary of hash_Sz basis: '
        !do i=1, Sz_subspace_size
        !    write(* ,*) i, hash_Sz(i)
        !end do

        !call H_XXX_subspace_fill_and_diag_fast(N_spin, J_spin, Sz_subspace_size, hash_Sz) !here only csr3

        !deallocate(hash, basis_rho_target, eigen_values, eigen_vectors, rho_reduced)
        !call cpu_time(finish_one_spin)
        !print *, "Finished for ", target_sz_spin
        !print *, " "
        !print '("Time Sz loop = ",f," seconds.")',finish_one_spin-start_one_spin

    !end do 
    
    deallocate (hash_Sz)

    !BELOW OLD VERSION OF THE CODE
    
    !creating a basis for N_spin system and whole basis vector -> 0000, 0001
    !print *, "Creating a basis for N_spin system and whole basis vector"
    !call H_create_basis_sz(N_spin, Sz_basis, basis_vector)
    !print *, "Basis created"
    !write(*,*) " "
    
    !preparing the file for the main loop
    !file_name1 = 'Entropy_results_' // trim(adjustl(N_spin_char)) // '_equal_div_feast_MACIEK_TEST.dat'
    !open (unit=1, file= trim(file_name1), recl=512)
    !write(1,*), "# Eigenvalue     Entropie    Spin_z    Sum_of_lambdas"
    
    !print*, "Entering main loop of the program"
                !spin_basis = dec2bin
    !size_of_sub_A = int(N_spin/2) !equally divided A and B systems
    !size_of_sub_B = int(N_spin - size_of_sub_A)
    !print *, "Size of subsystem A", size_of_sub_A
    !print *, "Size of subsystem B", size_of_sub_B
    !print *, " "

    !e_up = 20d0 
    !e_down = -20d0  
    

    !emin = e_down ! The lower ... &
    !emax =  e_up  !  and upper bounds of the interval to be searched for eigenvalues

    ! !$OMP PARALLEL DO
    !do k = 1, size(target_sz, 1) ! loop over all S_z basis
            !chosen target in this loop
            !call cpu_time(start_one_spin)
            !target_sz_spin = target_sz(k)

            ! here put loop for whole target_sz
            !call Hash_basis_with_target(N_spin, target_sz_spin , Sz_basis, basis_vector, hash, basis_rho_target)
            !print *, "Start entropy calculation for chosen target of Sz: ", target_sz_spin

            ! calulacting basis for subsystems A and B 

            !call H_XXX_block_diag_with_target_dense(N_spin, J_spin, hash)
            !print *, "LAPACK dense block matrix diagonalization for Sz: ", target_sz_spin
            !call H_XXX_block_csr_lapack_vec_diag_dense(N_spin, J_spin, no_of_nonzero, hash, eigen_values, eigen_vectors)

            !number_of_eigen = MAXVAl(hash)
            !print *, "Feast vectors filling for Sz: ", target_sz_spin
            !call H_XXX_block_feast_vec_fill(N_spin, J_spin, hash, no_of_nonzero)
            !print *, "Feast diagonalization for Sz: ", target_sz_spin
            ! subroutine H_XXX_block_feast_vec_diag(N_spin, e_up, e_down, number_of_eigen, no_of_nonzero, hash, e, x)
            !call H_XXX_block_feast_vec_diag(N_spin, e_up, e_down, number_of_eigen, no_of_nonzero, hash, eigen_values, eigen_vectors)

            !print *, "Loop of generating reduced density matrices and entropy calculation for Sz: ", target_sz_spin

            ! !$OMP PARALLEL DO
            !do i = 1, size(eigen_values, 1)
                !call Rho_reduced_calculation(N_spin, basis_rho_target, size_of_sub_A, size_of_sub_B, eigen_vectors, i, rho_reduced)
                !call Entropy_calculation(size_of_sub_A, size(rho_reduced, 1), rho_reduced, entropy_value, eigen_value_check)

                !print*, "This is entropy value for:"
                !print*, "Spin Z = ", target_sz_spin
                !print*, "Eigenvalue = ", eigen_values(i)
                !print*, "Entropy normalised = ", entropy_value
                !print*, "Lambdas check if 1.0 = ", eigen_value_check
                !write(1,*) eigen_values(i), ',' , entropy_value, ',',  target_sz_spin, ',', eigen_value_check

            !end do 
            ! !$OMP END PARALLEL DO

            !deallocate(hash, basis_rho_target, eigen_values, eigen_vectors, rho_reduced)
            !call cpu_time(finish_one_spin)
            !print *, "Finished for ", target_sz_spin
            !print *, " "
            !print '("Time Sz loop = ",f," seconds.")',finish_one_spin-start_one_spin
    ! end do 
    ! !$OMP END PARALLEL DO
    
    close(1)
    !deallocate(Sz_basis, basis_vector, target_sz, indices_Sz_basis_sorted)


    write(*,*) " "
    write(*,*) "Program executed with success"
    call cpu_time(finish)
    print '("Time = ",f," seconds.")',finish-start
    write(*,*) '------- END Heisenberg Program -------'

    

end program spin_code