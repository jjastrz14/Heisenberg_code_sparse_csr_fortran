module heisenberg 
    implicit none
    contains

    subroutine Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)
        use math_functions
        use timing_utilities
        implicit none
        integer, intent(in) :: N_spin
        double precision, intent(in) :: Sz_choice
        integer (8), intent(out) :: Sz_subspace_size
        integer, allocatable :: hash_Sz(:)

        integer :: N_spin_max, i, j, ind_2, Sz_choice_ind
        double precision :: Sz
        logical :: bool
        type(timer) :: calc_timer

        N_spin_max = 2**N_spin

        !from all N_spin_max we should be able to predict degeneracy of restricted S_z
        !1) we need to know Sz_choice_ind from Sz_choice and N_spin_max
        !Sz_choice = N_spin_max/2      -> Sz_choice_ind = 0
        !Sz_choice = N_spin_max/2-1    -> Sz_choice_ind = 1
        !Sz_choice = N_spin_max/2-2    -> Sz_choice_ind = 2
        write(*,*) 'choosen spin= ', Sz_choice

        Sz_choice_ind = int( N_spin/2.0d0 - Sz_choice )
        write(*,*) 'calculated integer Sz_choice_ind = ', Sz_choice_ind

        !testing binomialCoefficient(n, k, C)
        call binomialCoefficient(N_spin, Sz_choice_ind, Sz_subspace_size)

        !Sz_subspace_size = 4 ! standard factorial calcualated recursively/ logartihmically
        write(*,*) 'Sz_subspace_size = ', Sz_subspace_size

        !Sz_subspace_size = 1

        allocate(hash_Sz(Sz_subspace_size) )

        ! we can use e.g. btest function, which returns info if given bit is 0 (false) or 1 (true)
        !bool_T/F=btest(number, bit_number)
        ! 0 == false => Sz=Sz+1/2, 1= true> = Sz=Sz-1/2

        call calc_timer%start()

        !measure time hash Sz
        ind_2 = 1
        do i = 0 , N_spin_max-1
            !basis(i) = i !here we want Sz from this i
            !write(*,*) 'basis vector i =', i+1
            Sz = 0.0d0
            do j=0, N_spin-1
                bool = btest(i, j)
                !write(*,*) j, bool

                if (bool==.FALSE.) then
                    Sz=Sz+1.0d0/2.0d0
                elseif (bool==.TRUE.) then
                    Sz=Sz-1.0d0/2.0d0
                else
                    write(*,*) 'something wrong'
                    error stop
                end if
            end do

            !Sz_basis(i+1) = Sz
            if (Sz==Sz_choice) then
                hash_Sz(ind_2) = i+1
                ind_2 = ind_2 + 1
            end if

            !write(*,*) 'for basis vector i =', i, 'Sz = ', Sz
        end do

        write(*,*) 'internal hash_Sz size check', ind_2-1, 'vs', Sz_subspace_size
        write(*,*) 'remember that this vector can in principle be evaluated dynamically'

        call calc_timer%stop()
        write(*,*) 'Time of Sz subspace choice'
        call calc_timer%print_elapsed(time_unit%seconds, "seconds")
        call calc_timer%print_elapsed(time_unit%minutes, "minutes")
        call calc_timer%print_elapsed(time_unit%hours, "hours")

        !open (unit=40, file="hash_Sz_basis.dat", recl=512)
        !write(*,*) 'Summary of hash_Sz basis: '
        !do i=1,ind_2-1
        !    !write(* ,*) i, hash_Sz(i)
        !    write(40,*) i, hash_Sz(i)
        !end do
        !close(40)

    end subroutine Sz_subspace_choice


    subroutine H_XXX_filling(N_spin, J_spin, i, j, H_full_ij)
        implicit none
        integer, intent(in) :: N_spin, i, j
        double precision, intent(in) :: J_spin
        double precision, intent(out) :: H_full_ij
        integer :: N_spin_max  
        double precision :: Sp_site_ind_Sm_site_ind_p1_ij, Sm_site_ind_Sp_site_ind_p1_ij, Sz_site_ind_Sz_site_ind_p1_ij
        integer :: site_ind, m, n, p, q, r, s
        integer :: ind_i, ind_j, a_i, a_j, b_i, b_j, c_i, c_j
        double precision :: A_val, B_val, C_val

        ! N spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        N_spin_max = 2**N_spin

        H_full_ij = 0.0d0
                           
        do site_ind = 1, (N_spin-1) !here change for searching through the adjacency matrix for graphs
            m = 2**(site_ind-1)
            n = m
            p = 4 ! N_spin independent because this referes to S+S- operator for spin 1/2 
            q = 4
            r = 2**(N_spin-(site_ind+1))
            s = r
            !write(*,*) site_ind, m, n, p, q, r, s
        
            Sp_site_ind_Sm_site_ind_p1_ij = 0.0d0
            Sm_site_ind_Sp_site_ind_p1_ij = 0.0d0
            Sz_site_ind_Sz_site_ind_p1_ij = 0.0d0
            
            !do ind_i = 1, N_spin_max
            !do ind_j = 1, N_spin_max
            ind_i = i
            ind_j = j
                !write(*,*) 'site_ind, ind_i, ind_j = ',  site_ind, ind_i, ind_j
                !i = ind_i-1
                !j = ind_j-1
                a_i = FLOOR( (ind_i-1.0d0)/( (p*r)+0.0d0 ) )
                a_j = FLOOR( (ind_j-1.0d0)/( (q*s)+0.0d0 ) )
                b_i = MOD(FLOOR( (ind_i-1.0d0)/( r+0.0d0 ) ), p)
                b_j = MOD(FLOOR( (ind_j-1.0d0)/( s+0.0d0 ) ), q)
                c_i = MOD(ind_i-1, r)
                c_j = MOD(ind_j-1, s)
                                
                A_val = 0.0d0
                B_val = 0.0d0
                C_val = 0.0d0
            
                if (a_i == a_j) then
                    A_val = 1.0d0
                else
                    A_val = 0.0d0
                end if
            
                if (c_i == c_j) then
                    C_val = 1.0d0
                else
                    C_val = 0.0d0
                end if
            
                ! For S+S-:
                if (b_i == 1 .AND. b_j == 2) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if
                Sp_site_ind_Sm_site_ind_p1_ij = A_val*B_val*C_val
            
                ! For S-S+:
                if (b_i == 2 .AND. b_j == 1) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if                    
                Sm_site_ind_Sp_site_ind_p1_ij = A_val*B_val*C_val
            
                ! For SzSz:
                if (b_i == 0 .AND. b_j == 0) then
                    B_val = 0.250d0
                else if (b_i == 1 .AND. b_j == 1) then
                    B_val = -0.250d0
                else if (b_i == 2 .AND. b_j == 2) then
                    B_val = -0.250d0
                else if (b_i == 3 .AND. b_j == 3) then
                    B_val = 0.250d0
                else
                    B_val = 0.0d0
                end if

                Sz_site_ind_Sz_site_ind_p1_ij= A_val*B_val*C_val
            
            !end do
            !end do
        
            H_full_ij = H_full_ij + 1.0d0/2.0d0 *(Sp_site_ind_Sm_site_ind_p1_ij+Sm_site_ind_Sp_site_ind_p1_ij) + Sz_site_ind_Sz_site_ind_p1_ij
            
        end do


    end subroutine H_XXX_filling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!MPI REDISTRIBUTION OF HAMILTONIAN MATRIX FOR MULTI-NODE DIAGONALIZATION!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine Hamiltonian_diag_pfeast_multi_node_full_matrix(N_spin, J_spin, Sz_subspace_size, ia, ja, val_arr)
            use omp_lib
            use mkl_vsl
            use math_functions
            use timing_utilities
            implicit none
            include 'mpif.h'
    
            integer, intent(in) :: N_spin
            integer (8), intent(in) :: Sz_subspace_size
            double precision, intent(in) :: J_spin
            integer, allocatable, intent(inout) :: ia(:), ja(:)
            double precision, allocatable, intent(inout) :: val_arr(:)

            integer :: i, j, ind_i, ind_j, N_spin_max
            double precision :: H_full_ij, element_value_cutoff
            double precision :: norm 
            double precision, allocatable :: x(:,:)
            integer :: fpm(128)
            double precision :: emin, emax
            integer :: m0
            double precision :: epsout
            integer :: loop
            double precision, allocatable :: e(:)
            integer :: m
            double precision, allocatable :: res(:)
            CHARACTER*1 :: jobz, range, uplo
            integer :: il, iu, ldz, liwork, lwork, info, m_eig, n
            integer :: nL3, rank3, code
            integer, allocatable :: ia_local(:), ja_local(:)
            double precision, allocatable :: val_arr_local(:)
            double precision, allocatable :: gathered_eigenvalues(:), gathered_eigenvectors(:,:)
            integer :: rows_per_proc, start_row, end_row, local_nnz

            character(len=53) :: file_name   
            type(timer) :: calc_timer

            !this function was tested for the full matrix and should work (10 April 2025 - JJ)
    
            if (rank3 == 0) then

                call calc_timer%start()
                print *, "----------------------------------------"
                print *, "START FEAST diagonalisation"
                print *, "----------------------------------------"
                !debug: 
                do i=1,size(ia)
                    write(*,*) "index array: ", ia(i)
                end do

                do i=1, size(ja)
                    write(*,*) "ja and values array: ", ja(i), val_arr(i)
                end do 
            end if


            call MPI_COMM_SIZE(MPI_COMM_WORLD, nL3, code) !nL3 is the number of nodes!

            call pfeastinit(fpm, MPI_COMM_WORLD, nL3)  ! function specifying default parameters fpm of FEAST algorithm

            !MPI comunicator scheduling: 
            call MPI_COMM_RANK(fpm(49), rank3, code)
            
            rows_per_proc = Sz_subspace_size / nL3 !size of the problem divided by number of nodes

            if (rank3 == nL3 - 1) then
                rows_per_proc = Sz_subspace_size - (nL3 - 1) * rows_per_proc
            endif

            ! Calculate local row range for this process
            start_row = rank3 * (Sz_subspace_size / nL3) + 1
            end_row = start_row + rows_per_proc - 1

            ! Make sure we're not going past array bounds
            if (end_row > Sz_subspace_size) then
                end_row = Sz_subspace_size
            endif
            
            ! Calculate number of non-zeros in local portion
            local_nnz = ia(end_row + 1) - ia(start_row)

            ! Debug output to verify distribution
            write(*,*) 'Process ', rank3, ' handling rows ', start_row, ' to ', end_row
            write(*,*) 'Local non-zeros: ', local_nnz

            ! Allocate local arrays with verified sizes
            allocate(ia_local(rows_per_proc + 1), stat=info)  ! +1 for CSR format
            if (info /= 0) write(*,*) 'Error allocating ia_local on process ', rank3
            allocate(ja_local(local_nnz), stat=info)
            if (info /= 0) write(*,*) 'Error allocating ja_local on process ', rank3
            allocate(val_arr_local(local_nnz), stat=info)
            if (info /= 0) write(*,*) 'Error allocating val_arr_local on process ', rank3
            
            ! Fill local arrays
            ! Adjust ia indices to start from 1 in local numbering
            do i = 1, end_row - start_row + 2
                ia_local(i) = ia(start_row + i - 1) - ia(start_row) + 1
            end do
            
            ! Copy corresponding ja and values entries
            ja_local(1:local_nnz) = ja(ia(start_row):ia(end_row+1)-1)
            val_arr_local(1:local_nnz) = val_arr(ia(start_row):ia(end_row+1)-1)

            !int 64 has 64 bits
            !double precision has 64 bit
    
            !n=n     ! Sets the size of the problem
            !a=non_zero_array     ! Array containing the nonzero elements of the upper triangular part of the matrix A
        
            fpm(1) = 1
            fpm(2) = 12 !can be more, can be less
            fpm(3) = 12 !eps 10^-fpm(3)
            fpm(4) = 20 ! max number of feast loops
            fpm(5) = 0 !initial subspace
            fpm(6) = 0! stopping criterion
            fpm(7) = 5 !Error trace sigle prec stop crit
            fpm(14) = 0 ! standard use of feast
            fpm(27) = 1 !check input matrices
            fpm(28) = 1 !check if B is positive definite?
            
            uplo='F' !'U' ! If uplo = 'U', a stores the upper triangular parts of A.
            n = Sz_subspace_size ! Sets the size of the problem
            !Intervals for 10 lowest eigenstates for 1D chain of H_XXX with NN hoping
            emin = -20!-0.4465d0 * N_spin + 0.1801d0
            !emax = -0.49773d0 * N_spin + 2.10035d0 !it might be more adjusted to the possible number of eigenvalues to be found
            emax = 20!-0.49773d0 * N_spin + 2.10030d0
            
            if (rank3 == 0) then
                write(*,*) "Calculated lower bound: ", emin
                write(*,*) "Calculated upper bound: ", emax
            end if 

            !m0 = 20 !Sz_subspace_size !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n.
            m0 = Sz_subspace_size
            !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax].
            !If the initial guess is wrong, Extended Eigensolver routines return info=3.

            allocate( x(rows_per_proc,m0), e(m0), res(m0) )
            write(*,*) 'Memory study 4: n,m0: ', n, m0
            !write(*,*) 'Memory study 4: x: ', kind(x(1,1)) * size(x) , 'bytes',  kind(x(1,1)) * size(x)/1024.0/1024.0, 'MB'
            !write(*,*) 'Memory study 4: e: ', kind(e(1)) * size(e) , 'bytes',  kind(e(1)) * size(e)/1024.0/1024.0, 'MB'
            !write(*,*) 'Memory study 4: res: ', kind(res(1)) * size(res) , 'bytes',  kind(res(1)) * size(res)/1024.0/1024.0, 'MB'

            ! Debug output before FEAST
            write(*,*) 'Process ', rank3, ' ready for FEAST with:'
            write(*,*) '  Local rows: ', rows_per_proc
            write(*,*) '  Local non-zeros: ', local_nnz
            !write(*,*) '  Memory allocated for x:', rows_per_proc * m0 * 8, ' bytes'
            do i=1,size(ia_local)
                write(*,*) "index array: ", ia_local(i)
            enddo

            do i=1, size(ja_local)
                write(*,*) "ja and values array: ", ja_local(i), val_arr_local(i)
            enddo 

            call MPI_BARRIER(MPI_COMM_WORLD, code)

            write(*,*) 'Before dfeast_scsrev... '
            call pdfeast_scsrev(uplo, rows_per_proc, val_arr_local, ia_local, ja_local, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
            write(*,*) 'eps_out= ', epsout
            write(*,*) 'loop= ', loop
            write(*,*) ' dfeast_scsrev info=', info
            write(*,*) 'After  dfeast_scsrev... '
    
            if (info /= 0) then
                write(*,*) 'problem with  dfeast_scsrev, info=', info
            end if
    
            write(*,*) 'Process ', rank3, 'pdfeast_scsrev eigenvalues found= ', m
            
            ! console print of eigenvals for debug
            
            write(*,*) 'Process ', rank3, ' eigenvalue:'
            do i=1,m
                norm = 0.0
                do j=1, n
                    norm = norm + x(j,i)*(x(j,i))
                end do
                write(*,*) 'Eigenvalue', i
                write(*,*) i, e(i), norm, "node: ", rank3
            end do
            write(*,*) 'Process ', rank3, ' eigenvalues printed'

            call MPI_BARRIER(MPI_COMM_WORLD, code)

            if (rank3 == 0) then
                allocate(gathered_eigenvalues(m*nL3))
                allocate(gathered_eigenvectors(m*nL3,n))
                write(*,*) "allocation at master succesful!"
            end if

            !to ponizej nie potrzebne, bo kazdy procek ma swoje eigenvalues, nie wiem jak ogarnac eigenvecs
            ! Gather all eigenvalues and eigenvectors to rank3 == 0
            call MPI_GATHER(e, m, MPI_DOUBLE_PRECISION, &
            gathered_eigenvalues, m, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, code)

            call MPI_GATHER(x, m*n, MPI_DOUBLE_PRECISION, &
            gathered_eigenvectors, m*n, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, code)
  
            deallocate(ia_local, ja_local, val_arr_local)

            deallocate(x, e, res)

            write(*,*) "After MPI gather"

            call MPI_BARRIER(MPI_COMM_WORLD, code)

            if (rank3 == 0) then !only rank3=0 needs to do it

                write(file_name, '(A,I0,A)') 'Eigenvalues_results_', N_spin, '_feast_test_full.dat'
                open(10, file=trim(file_name), recl = 512)
                
                !save to file done by single node
                do i=1,m*nL3
                    norm = 0.0
                    do j=1, n
                       norm = norm + gathered_eigenvectors(j,i)*(gathered_eigenvectors(j,i))
                    end do
                    write(10,*) i, gathered_eigenvalues(i), norm
                end do

                call calc_timer%stop()
                print *, "----------------------------------------"
                print *, "END FEAST diagonalisation:"
                call calc_timer%print_elapsed(time_unit%seconds, "seconds")
                call calc_timer%print_elapsed(time_unit%minutes, "minutes")
                call calc_timer%print_elapsed(time_unit%hours, "hours")
                print *, "----------------------------------------"
                call calc_timer%reset()
                close(10)
            endif    
            
        end subroutine Hamiltonian_diag_pfeast_multi_node_full_matrix


        subroutine H_XXX_MPI_OpenMP_pfeast(hash_Sz, Sz_subspace_size,  N_spin, J_spin)        
            use omp_lib    
            implicit none    
            include 'mpif.h'
            integer(8), intent(in) :: Sz_subspace_size, N_spin
            double precision, intent(in) :: J_spin
            integer(8), allocatable :: hash_Sz(:)
            integer(8) :: N, MPI, perfect_N, elements_left, ind, ind2, mpi_row_start, ind_mpi, rows_per_mpi, el_in_row, el_to_be_burn, total_rows_burn, sum_check, no_of_els_to_calculate, row_start, row_end, id_ind, window_size, windows_per_row, row, ind_inside_window, ind_i, ind_j, window, last_index, ia_temp_ind, temp_ind, no_of_ja_val_el, i
            integer(8), allocatable :: mpi_rows_collected(:,:), ia_temp_vectors(:), ja_temp_vectors(:)
            integer :: ierror, size_of_cluster, local_process_ID
            double precision, allocatable :: val_temp_vectors(:)
            double precision :: H_full_ij,  element_value_cutoff
            character(len=32) :: loc_id_info
            
            open (unit=10, file="mpi_rows_collected_test_24Mar2025.dat", recl=512)
            !open (unit=11, file="mpi_H_31Mar2025.dat", recl=512)
            !open (unit=12, file="mpi_Hels_31Mar2025.dat", recl=512)
            
            element_value_cutoff =  10.0d0**(-8.00d0) ! only larger elements are taken as non-zero
                    
            call MPI_COMM_SIZE(MPI_COMM_WORLD, size_of_cluster, ierror) !Returns the size of the group associated with a communicator.
            if (ierror .NE. 0) write(*,*) 'error after MPI_COMM_SIZE ', ierror
            write(*,*) 'size of the cluster is = ', size_of_cluster
        
            call MPI_COMM_RANK(MPI_COMM_WORLD, local_process_ID, ierror) !Determines the rank of the calling process in the communicator.
            if (ierror .NE. 0) write(*,*) 'error after MPI_COMM_RANK ', ierror
            write(*,*) 'local_process_ID (rank) is = ', local_process_ID

            write(*,*) "List of processes before barrier", local_process_ID
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            write(*,*) "List of processes after barrier", local_process_ID
            
            !test delete later
            do i=1,size_of_cluster
            if (local_process_ID == i-1) then
                write(*,*) 'local_process_ID = ', local_process_ID, 'hash_Sz(1-3) = ', hash_Sz(1), hash_Sz(2), hash_Sz(3)
            end if
            end do
                
            ! so now actual code for ia, ja, val_arr save on disk (then we can try to re-allocate)
                
            !no of rows per mpi process
            !last number of rows might be different, eg. 20/3 gives 7 7 6(!)
            !N_mat x N_mat = Sz_subspace_size x Sz_subspace_size
            rows_per_mpi = NINT((Sz_subspace_size+0.0d0)/(size_of_cluster+0.0d0))
            write(*,*) 'nearest int of Sz_subspace_size/size_of_cluster' ,local_process_ID, rows_per_mpi
                
            if (local_process_ID == size_of_cluster - 1) then
                     rows_per_mpi = Sz_subspace_size - (size_of_cluster - 1) * rows_per_mpi
            endif
                
            write(*,*) 'test of rows_per_proc', local_process_ID, rows_per_mpi
                
                
            ! row_start and row_end per mpi process
                
            ! Calculate local row range for this process
            start_row = 1 +  local_process_ID * NINT((Sz_subspace_size+0.0d0)/(size_of_cluster+0.0d0))
            end_row = start_row + rows_per_mpi - 1
            !    
            ! Make sure we're not going past array bounds
            if (end_row .GT. Sz_subspace_size) then
                    error stop
            endif
            
            size_of_cluster = 3

            N=Sz_subspace_size !84756;
            write(*,*) 'N=Sz_subspace_size', N
            MPI = size_of_cluster
            allocate(mpi_rows_collected(MPI,2))
            !%perfect_N = floor((N/2*(N+1))/MPI);
            perfect_N = NINT((N/2.0d0*(N+1.0d0))/(MPI + 0.0d0));
            write(*,*) 'perfect_N=', perfect_N
        
            !below calculation how many rows per mpi process
            elements_left = perfect_N
            ind = 0
            mpi_row_start = 0
            
            do ind_mpi=1, (MPI-1)
                rows_per_mpi = 1
                loop1: do ind=mpi_row_start, (N-1)        
                    el_in_row = N - ind
                    el_to_be_burn = elements_left - el_in_row
                    if (el_to_be_burn .GT. (el_in_row-1)) then
                        rows_per_mpi=rows_per_mpi+1
                        elements_left = el_to_be_burn
                    else            
                        exit loop1
                    end if
                end do loop1
                mpi_row_start = rows_per_mpi    
                el_to_be_burn=el_to_be_burn+perfect_N
                elements_left = el_to_be_burn
            
                mpi_rows_collected(ind_mpi,1)=rows_per_mpi
                write(*,*) 'rows_per_mpi = ', rows_per_mpi
            end do
            
            ! last process must calcualte all remaining elements
            total_rows_burn = 0
            do ind_mpi=1, (MPI-1 )
                total_rows_burn = total_rows_burn+mpi_rows_collected(ind_mpi,1)
            end do
            mpi_rows_collected(MPI,1)=N-total_rows_burn
            write(*,*) 'rows_per_mpi last= ', N-total_rows_burn
            
            ! first column - number of rows, second column = number of elements
            !workload balance test
            ind=0
            do ind_mpi=1,MPI
                no_of_els_to_calculate = 0
                do ind2=1, mpi_rows_collected(ind_mpi,1)
                    no_of_els_to_calculate = no_of_els_to_calculate+(N-ind)
                    ind=ind+1
                end do
                mpi_rows_collected(ind_mpi,2) = no_of_els_to_calculate
            end do
            
            sum_check = 0
            do ind_mpi=1,MPI
                sum_check = sum_check + mpi_rows_collected(ind_mpi,2)
                write(*,*) 'ind_mpi. no_of_el' , ind_mpi, mpi_rows_collected(ind_mpi,2)
            end do
            
            write(*,*) 'sum_check .EQ. (N/2*(N+1))', sum_check, (N*(N+1))/2
            
            if (  sum_check .EQ. (N*(N+1))/2 ) then
                write(*,*) 'rows checked, seems ok'
            else
                write(*,*) 'somethings wrong'
                error stop
            end if
            
            !write to file
            do ind = 1, mpi        
                write(10,*) mpi_rows_collected(ind,1), mpi_rows_collected(ind,2)        
            end do


            !!!!!!!!      CREATION OF IA JA AND VALUES ARRAY               !!!!!!!!!!!!!
            
            ! after procedure above we have access to mpi_rows_collected(ind_mpi,1-2)
            
            !local_process_ID = 0 -> 3 rows (1,2,3)
            !local_process_ID = 1 -> 5 rows (4,5,6,7,8)
            !local_process_ID = 2 -> 12 rows (9,10 - 20) = (3+5+1, 3+5+12rows)
            
            !write(*,*) 'test if 9', CEILING(8.3)

            !call once again the rank of the process we are in
            call MPI_COMM_RANK(MPI_COMM_WORLD, local_process_ID, ierror) !Determines the rank of the calling process in the communicator.
            if (ierror .NE. 0) write(*,*) 'error after MPI_COMM_RANK ', ierror
            write(*,*) 'local_process_ID (rank) is = ', local_process_ID
            
            local_process_ID = 0
            
            write (loc_id_info,'(i0)') local_process_ID
            loc_id_info = adjustl(loc_id_info)
            open (unit=15, file='ia_local_'//trim(loc_id_info)//'.dat', recl=512)
            open (unit=16, file='ja_local_'//trim(loc_id_info)//'.dat', recl=512)
            open (unit=17, file='val_local_'//trim(loc_id_info)//'.dat', recl=512)
            
            window_size = 8 ! should be omp_procs * x
            
            !$OMP PARALLEL PRIVATE(omp_id)
            omp_id = omp_get_thread_num()
            write(*,*) 'OpenMP thread check', omp_id
            !$OMP END PARALLEL

            !windows_size = omp_id
            
            if(window_size.GE.Sz_subspace_size) then
                write(*,*) 'too large window size'
                error stop
            end if
            
            !allocation of 3 vectors matrix
            allocate (ia_temp_vectors(mpi_rows_collected(local_process_ID + 1,1)+1), ja_temp_vectors(window_size), val_temp_vectors(window_size) )
            
            row_start = 1
            do ind = 1, local_process_ID
                row_start = row_start + mpi_rows_collected(ind ,1)
            end do
            row_end = row_start + mpi_rows_collected(local_process_ID + 1,1) -1 !<- how many rows per mpi process
            
            !write(*,*) 'row_start, row_end', row_start, row_end looks ok 28Mar2025
            
            !ind = 1
            write(*,*) 'row_start, row_end', row_start, row_end
            no_of_ja_val_el = 0
            ia_temp_ind = 0
            temp_ind =1
            do row = row_start, row_end
                ia_temp_vectors(temp_ind) = ia_temp_ind+1
                temp_ind = temp_ind + 1
                windows_per_row = CEILING( (Sz_subspace_size-row+1.0d0)/(window_size+0.0d0) )
                write(*,*) 'row, windows_per_row', row, windows_per_row
                do window = 1, windows_per_row - 1 !becasue last one will be special
                    ja_temp_vectors  = 0
                    val_temp_vectors = 0.0d0
                    
                    !$OMP PARALLEL DO PRIVATE(ind_inside_window, ind_i, ind_j, H_full_ij) SHARED(window_size, row, hash_Sz, window, N_spin, J_spin,element_value_cutoff, ja_temp_vectors, val_temp_vectors, ia_temp_ind) default(none)
                    do ind_inside_window = 1, window_size
                        ind_i = hash_Sz(row)
                        ind_j = hash_Sz( (row-1) + ind_inside_window + (window-1)*window_size )
                        call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)

                        if ( (abs(H_full_ij) .GT. element_value_cutoff) .AND. (ind_i .NE. ind_j) ) then
                            !collect non-zero in temp_vectors
                            ja_temp_vectors(ind_inside_window) = (row-1) + ind_inside_window + (window-1)*window_size
                            val_temp_vectors(ind_inside_window) = H_full_ij
                        else if (ind_i == ind_j) then
                            ja_temp_vectors(ind_inside_window) = (row-1) + ind_inside_window + (window-1)*window_size
                            val_temp_vectors(ind_inside_window) = H_full_ij
                        end if
                    
                    end do
                    !$OMP END PARALLEL DO
                    
                    ! since ja_temp_vectors, val_temp_vectors will be zeroed in a moment, we should save/collect them somehow here
                    
                    do i=1,window_size
                        if ( (abs(val_temp_vectors(i)) .GT. element_value_cutoff)) then
                            write(16,*) ja_temp_vectors(i)
                            write(17,*) val_temp_vectors(i)
                            ia_temp_ind = ia_temp_ind + 1
                            no_of_ja_val_el = no_of_ja_val_el + 1
                        end if
                    end do
                end do
                
                !last window
                window = windows_per_row
                    ja_temp_vectors  = 0
                    val_temp_vectors = 0.0d0
                    
                    if (modulo((Sz_subspace_size-row+1), window_size) == 0) then
                        last_index = window_size
                    else
                        last_index = modulo((Sz_subspace_size-row+1), window_size)
                    end if
                    
                !$OMP PARALLEL DO PRIVATE(ind_inside_window, ind_i, ind_j, H_full_ij) SHARED(window_size, last_index, row, hash_Sz, window, N_spin, J_spin,element_value_cutoff, ja_temp_vectors, val_temp_vectors, ia_temp_ind) default(none)
                do ind_inside_window = 1, last_index !(in last window set to mod(x,y) )
                    !write(11,*) row, window, ind_inside_window
                    ind_i = hash_Sz(row)
                    ind_j = hash_Sz( (row-1) + ind_inside_window + (window-1)*window_size )
                    !write(*,*) 'i,j', ind_i, ind_j
                    call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
                    
                    if ( (abs(H_full_ij) .GT. element_value_cutoff) .AND. (ind_i .NE. ind_j) ) then
                        !collect non-zero in temp_vectors
                        ja_temp_vectors(ind_inside_window) =  (row-1)+ (window-1)*window_size + ind_inside_window
                        val_temp_vectors(ind_inside_window) = H_full_ij

                    else if (ind_i == ind_j) then
                        ja_temp_vectors(ind_inside_window) = (row-1) + ind_inside_window + (window-1)*window_size
                        val_temp_vectors(ind_inside_window) = H_full_ij

                    end if
                    
                end do
                !$OMP END PARALLEL DO
                
                ! since ja_temp_vectors, val_temp_vectors will be zeroed in a moment, we should save/collect them somehow here
                do i=1, last_index
                    if ( (abs(val_temp_vectors(i)) .GT. element_value_cutoff)) then
                        write(16,*) ja_temp_vectors(i)
                        write(17,*) val_temp_vectors(i)
                        no_of_ja_val_el = no_of_ja_val_el + 1
                        ia_temp_ind = ia_temp_ind + 1
                    end if
                end do
            
            end do
            
            ia_temp_vectors(temp_ind) = ia_temp_ind+1
            
            !here we must collect temp_vec to final ia, ja, val_arr realloc vs save on disk
            
            do i = 1, mpi_rows_collected(local_process_ID + 1,1) + 1
                write(*,*) 'i, ia_local(i)', i, ia_temp_vectors(i)
                ! this should be saved to ia_local file
                write(15,*) ia_temp_vectors(i)
            end do
        
            
            deallocate (ia_temp_vectors, ja_temp_vectors, val_temp_vectors )
            deallocate (mpi_rows_collected)
            close(10)
            !close(11)
            !close(12)
            
            close(15)
            close(16)
            close(17)
            
            end subroutine H_XXX_MPI_OpenMP_pfeast

end module heisenberg 