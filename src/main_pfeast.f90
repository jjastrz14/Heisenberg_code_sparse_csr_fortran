program spin_code
    use heisenberg
    use tests_module
    use timing_utilities
    use mpi
    implicit none

    integer (8) :: N_spin, N_spin_max, Sz_subspace_size, size_ia, size_ja, size_val, number_of_rows_local_MPI
    double precision :: J_spin, Sz_choice
    character(len = 12) :: N_spin_char, J_spin_char
    integer (8), allocatable :: hash_Sz(:) !check if it should be integer(8)
    integer :: local_process_ID, size_of_cluster, code, ierror
    type(timer) :: calc_timer_main

    ! Initialize the MPI environment
    call MPI_INIT(ierror)
    if (ierror .NE. 0) write(*,*) 'error after MPI_INIT ', ierror
    call MPI_COMM_RANK(MPI_COMM_WORLD, local_process_ID, ierror)
    if (ierror .NE. 0) write(*,*) 'error after first MPI_COMM_RANK ', ierror
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size_of_cluster, ierror)
    if (ierror .NE. 0) write(*,*) 'error after first MPI_COMM_SIZE ', ierror

    ! Process command line arguments
    If(command_argument_count().NE.2) Then
        write(*,*)'Error, Only N (integer) and J (double precision) is required'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierror)  ! Terminate all processes if error
    endif 

    ! Get and parse command line arguments
    call get_command_argument(1, N_spin_char)
    call get_command_argument(2, J_spin_char)
    read(N_spin_char, *) N_spin
    read(J_spin_char, *) J_spin

    !test for MPI+OpenMP
    call MPI_plus_OpenMP_test()

    !Sz_subspace calucation and hash_Sz generation
    call Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)

    ! Only rank 0 performs the initial setup and calculations
    if (local_process_ID == 0) then
        call calc_timer_main%start()

        ! Initial output - only from rank 0
        write(*,*) '------- START Heisenberg Program -------'
        write(*,*) 'Calculation of Heisenberg chain for N =', N_spin, 'and J = ', J_spin 
        write(*,*) ' '
    endif

    ! Calculate problem size and parameters
    N_spin_max = 2**N_spin 
    
    ! Choose subspace based on even/odd N_spin
    if (mod(N_spin,2) == 0.0d0) then 
        Sz_choice = 0.0d0 
    else 
        Sz_choice = 0.5d0  
    endif 
    
    ! These build the Hamiltonian matrix in CSR format
    call Hamiltonian_MPI_OpenMP_pfeast(hash_Sz, Sz_subspace_size,  N_spin, J_spin, size_ia, size_ja, size_val, number_of_rows_local_MPI) 

    ! Subroutine to PFEAST diagonalization
    call Hamiltonian_pfeast_diagonalisation(N_spin, Sz_subspace_size, size_ia, size_ja, size_val, number_of_rows_local_MPI)

    write(*,*) 'List of processes before MPI barrier', local_process_ID
    call MPI_BARRIER(MPI_COMM_WORLD, ierror) !synchronize all processes
    if (ierror .NE. 0) write(*,*) 'error after MPI_BARRIER ', ierror
    write(*,*) 'List of processes after MPI barrier', local_process_ID

    if (local_process_ID == 0) then
        write(*,*) " "
        write(*,*) "Program executed with success"
        call calc_timer_main%stop()
        call calc_timer_main%print_elapsed(time_unit%seconds, "seconds")
        call calc_timer_main%print_elapsed(time_unit%minutes, "minutes")
        call calc_timer_main%print_elapsed(time_unit%hours, "hours")
    endif
    
    call MPI_FINALIZE(ierror)
    if (ierror .NE. 0) write(*,*) 'error after MPI_FINALIZE ', ierror
    write(*,*) '------- END Heisenberg Program -------'

end program spin_code