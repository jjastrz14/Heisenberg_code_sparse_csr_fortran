program spin_code
    use heisenberg
    use tests_module
    use timing_utilities
    implicit none

    integer :: N_spin, N_spin_max 
    integer (8) :: Sz_subspace_size
    double precision :: J_spin, Sz_choice
    character(len = 12) :: N_spin_char, J_spin_char
    integer, allocatable :: hash_Sz(:), ia(:), ja(:)
    double precision, allocatable :: target_sz(:), val_arr(:)
    type(timer) :: calc_timer_main

    call calc_timer_main%start()
    If(command_argument_count().NE.2) Then
        write(*,*)'Error, Only N (integer) and J (double precision) is required, program stopped'
        stop 
    endif 

    call get_command_argument(1, N_spin_char)
    call get_command_argument(2, J_spin_char)

    read(N_spin_char, *) N_spin
    read(J_spin_char, *) J_spin

    !from test_module
    write(*,*) '------- START test modules -------'

    ! call mmm_csr_test()
    !call omp_mkl_small_test()
    ! call test_permutation_H_for_4_sites()
    ! call sparse_zfeast_test()
    !call sparse_dfeast_test()

    write(*,*) '------- END test modules -------'

    !N_spin = 4
    !J_spin = 1 

    write(*,*) '------- START Heisenberg Program -------'
    write(*,*) 'Calculation of Heisenberg chain for N =', N_spin , 'and J = ', J_spin 
    write(*,*) ' '

    N_spin_max = 2**N_spin 
    
    !Sz_choice = 0.0d0 !integer counted from Sz_max (in a sense that Sz_choice = 1 means eg for N_spin=4, Sz_max=2 Sz=2)

    !Choose always biggest subspace of the Hamiltonian
    if (mod(N_spin,2) == 0.0d0) then !0 subspace for even 0.5 supspace for odd
        Sz_choice = 0.0d0 
    else 
        Sz_choice = 0.5d0  
    endif 
    
    call Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)

    call Hamiltonian_fill_open_mp(N_spin, J_spin, Sz_subspace_size, hash_Sz, ia, ja, val_arr)

    call Hamiltonian_diag_feast(N_spin, J_spin, Sz_subspace_size, ia, ja, val_arr)

    write(*,*) " "
    write(*,*) "Program executed with success"
    call calc_timer_main%stop()
    call calc_timer_main%print_elapsed(time_unit%seconds, "seconds")
    call calc_timer_main%print_elapsed(time_unit%minutes, "minutes")
    call calc_timer_main%print_elapsed(time_unit%hours, "hours")
    write(*,*) '------- END Heisenberg Program -------'

    

end program spin_code