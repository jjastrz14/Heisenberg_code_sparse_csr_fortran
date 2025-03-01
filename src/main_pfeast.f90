program spin_code
    use heisenberg
    use tests_module
    use timing_utilities
    use mpi
    implicit none

    integer :: N_spin, N_spin_max, code, ierror
    integer (8) :: Sz_subspace_size
    double precision :: J_spin, Sz_choice
    character(len = 12) :: N_spin_char, J_spin_char
    integer, allocatable :: ia(:), ja(:)
    integer, allocatable :: hash_Sz(:) !check if it should be integer(8)
    double precision, allocatable :: target_sz(:), val_arr(:)
    integer :: rank, nprocs
    
    integer :: ja_size, val_size
    type(timer) :: calc_timer_main

    ! Initialize the MPI environment
    call MPI_INIT(code)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, code)

    !test for MPI+OpenMP
    call MPI_plus_OpenMP_test()

    !Sz_subspace calucation and hash_Sz generation
    call Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)
    call Hamiltonian_mpi_redistribution(N_spin, J_spin, Sz_subspace_size, hash_Sz)

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    write(*,*) 'EVRYTHING BELOW TO DEBUG'
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    ! Only rank 0 performs the initial setup and calculations
    if (rank == 0) then
        call calc_timer_main%start()
        
        ! Process command line arguments
        If(command_argument_count().NE.2) Then
            write(*,*)'Error, Only N (integer) and J (double precision) is required'
            call MPI_ABORT(MPI_COMM_WORLD, 1, code)  ! Terminate all processes if error
        endif 

        ! Get and parse command line arguments
        call get_command_argument(1, N_spin_char)
        call get_command_argument(2, J_spin_char)
        read(N_spin_char, *) N_spin
        read(J_spin_char, *) J_spin

        ! Initial output - only from rank 0
        write(*,*) '------- START Heisenberg Program -------'
        write(*,*) 'Calculation of Heisenberg chain for N =', N_spin, 'and J = ', J_spin 
        write(*,*) ' '

        ! Calculate problem size and parameters
        N_spin_max = 2**N_spin 
        
        ! Choose subspace based on even/odd N_spin
        if (mod(N_spin,2) == 0.0d0) then 
            Sz_choice = 0.0d0 
        else 
            Sz_choice = 0.5d0  
        endif 
        
        ! Sequential calculations - only on rank 0
        ! These build the Hamiltonian matrix in CSR format
        call Hamiltonian_fill_open_mp(N_spin, J_spin, Sz_subspace_size, hash_Sz, ia, ja, val_arr)
        
        !Sz_subspace_size = 4
        !allocate(ia(5), ja(14), val_arr(14))
        !4x4 matrix, remember to adjust the number of nodes to the size of the problem
        !ia = (/1, 4, 8, 12, 15/)
        !ja = (/1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4/)
        !val_arr = (/2.0d0, -1.0d0, -1.0d0, -1.0d0, 3.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0, 3.0d0, -1.0d0, -1.0d0, -1.0d0, 2.0d0/)

        ! Calculate sizes for MPI communication
        write(*,*) "ia_size", size(ia)
        ja_size = ia(size(ia)) - 1  ! Size of ja array in CSR format
        write(*,*) "ja_size", ja_size
        val_size = ja_size      ! val_arr has same size as ja
    endif

    ! Broadcast essential parameters to all processes
    call MPI_BCAST(N_spin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
    call MPI_BCAST(J_spin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, code)
    call MPI_BCAST(Sz_subspace_size, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, code)
    call MPI_BCAST(ja_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
    call MPI_BCAST(val_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)

    ! Non-root processes allocate their arrays
    if (rank /= 0) then
        allocate(ia(Sz_subspace_size + 1))
        allocate(ja(ja_size))
        allocate(val_arr(val_size))
    endif

    ! Broadcast CSR format arrays to all processes
    call MPI_BCAST(ia, Sz_subspace_size + 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
    call MPI_BCAST(ja, ja_size, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
    call MPI_BCAST(val_arr, val_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, code)

    ! Now all processes have the matrix data and can participate in parallel diagonalization
    call Hamiltonian_diag_pfeast_multi_node_upper_train_matrix(N_spin, J_spin, Sz_subspace_size, ia, ja, val_arr)

    if (rank == 0) then
        write(*,*) " "
        write(*,*) "Program executed with success"
        call calc_timer_main%stop()
        call calc_timer_main%print_elapsed(time_unit%seconds, "seconds")
        call calc_timer_main%print_elapsed(time_unit%minutes, "minutes")
        call calc_timer_main%print_elapsed(time_unit%hours, "hours")
        write(*,*) '------- END Heisenberg Program -------'
    endif

    if (allocated(ia)) deallocate(ia)
    if (allocated(ja)) deallocate(ja)
    if (allocated(val_arr)) deallocate(val_arr)
    if (rank == 0 .and. allocated(hash_Sz)) deallocate(hash_Sz)

    call MPI_FINALIZE(code)

end program spin_code