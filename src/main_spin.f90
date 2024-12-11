module heisenberg 
    implicit none
    contains

    subroutine Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)
        use math_functions
        implicit none
        integer, intent(in) :: N_spin
        double precision, intent(in) :: Sz_choice
        integer (8), intent(out) :: Sz_subspace_size
        integer, allocatable :: hash_Sz(:)

        integer :: N_spin_max, i, j, ind_2, Sz_choice_ind
        double precision :: Sz
        logical :: bool

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


    subroutine Hamiltonian_fill_diag_open_mp(N_spin, J_spin, Sz_subspace_size, hash_Sz)
        use omp_lib
        use mkl_vsl
        use math_functions
        use timing_utilities
        implicit none

        integer, intent(in) :: N_spin
        integer (8), intent(in) :: Sz_subspace_size
        double precision, intent(in) :: J_spin
        integer, allocatable :: hash_Sz(:), list_of_ind_2(:,:), ia(:), ja(:), open_mp_counter(:)
        !integer(2), allocatable :: list_of_ind(:,:) ! Signed integer value from -32,768 to 32,767
        logical(1), allocatable :: list_of_ind_bool(:)
        integer :: i, j, ind_i, ind_j, N_spin_max,  size_of_1D_list_for_sweep, &
            ja_val_arr_size, ind_temp, ind_temp_2, omp_id, threads_max !ind_Sz_1, ind_Sz_2
        integer(8), allocatable :: list_of_ind(:,:) ! Signed integer value from -32,768 to 32,767
        integer(8) :: ind_3, size_of_list, ind_Sz_1, ind_Sz_2 !8 byte = 64 bit Signed integer value from -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
        double precision, allocatable :: val_arr(:)
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
        integer(8) :: start_count, end_count, count_rate
        double precision :: elapsed_time 
        character(len = 12) :: N_spin_char
        character(len=53) :: file_name   
        type(timer) :: calc_timer

        write(file_name, '(A,I0,A)') 'Eigenvalues_results_', N_spin, '_feast.dat'
        open(10, file=trim(file_name), access='stream')

        print *, "----------------------------------------"
        print *, "START CSR filling"
        print *, "----------------------------------------"
        
        call calc_timer%start()
        write(*,*) 'Memory study 0: hash_Sz: ', kind(hash_Sz(1)) * size(hash_Sz) , 'bytes', kind(hash_Sz(1)) * size(hash_Sz)/1024.0/1024.0, 'MB'

        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        if (N_spin .GT. 63) error stop !integer overflow predicted then
        N_spin_max = 2**N_spin
        !if(N_spin_max .GT. 9223372036854775807) error stop ! definite integer overflow
        !write(*,*) 'be carefull about N_spin_max max size - large integer possibility'
        write(*,*) 'test for N_spin_max', N_spin_max

        ! 1D list for fast sweep
        size_of_list = Sz_subspace_size*(Sz_subspace_size + 1  ) / 2
        allocate( list_of_ind_bool(size_of_list) ) ! we perhaps could re-think excistence of this list
        list_of_ind_bool = .FALSE.
        
        write(*,*) 'Memory study 1: size_of_list: ', size_of_list
        write(*,*) 'Memory study 1: list_of_ind_bool: ', kind(list_of_ind_bool(1)) * size(list_of_ind_bool) , 'bytes', kind(list_of_ind_bool(1)) * size(list_of_ind_bool)/1024.0/1024.0, 'MB'
        
        ! fast 1D sweep calculating size of ja and val_arr vectors
        write(*,*) 'OpenMP processors check'
        threads_max = omp_get_max_threads()
        write(*,*) 'OpenMP thread max num check', threads_max

        allocate(open_mp_counter(threads_max) )
        open_mp_counter = 0

        omp_id = omp_get_thread_num()
        !open_mp_counter(omp_id+1) = open_mp_counter(omp_id+1) + 1
        !write(*,*) 'OpenMP thread check', omp_id

        element_value_cutoff =  10.0d0**(-8.00d0) ! only larger elements are taken as non-zero
        !ja_val_arr_size = 0
     
        !Open mp loop 
        !$OMP PARALLEL DO PRIVATE(ind_3, omp_id, ind_Sz_1, ind_Sz_2, ind_i, ind_j, H_full_ij) SHARED(size_of_list, open_mp_counter, list_of_ind_bool, hash_Sz, N_spin, J_spin, element_value_cutoff, Sz_subspace_size) default(none)
        do ind_3 = 1, size_of_list
            omp_id = omp_get_thread_num()
            call index_map(Sz_subspace_size, ind_3, ind_Sz_1, ind_Sz_2)
            !ind_Sz_1 = list_of_ind(ind_3,1)
            !ind_Sz_2 = list_of_ind(ind_3,2)
            if (ind_Sz_1 == ind_Sz_2) then
                !ja_val_arr_size = ja_val_arr_size + 1
                open_mp_counter(omp_id+1) = open_mp_counter(omp_id+1) + 1
                list_of_ind_bool(ind_3) = .TRUE.
            else if (ind_Sz_2 .GT. ind_Sz_1) then
                ind_i = hash_Sz(ind_Sz_1)
                ind_j = hash_Sz(ind_Sz_2)
                call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
                if ( abs(H_full_ij) .GT. element_value_cutoff ) then
                    !ja_val_arr_size = ja_val_arr_size + 1
                    open_mp_counter(omp_id+1) = open_mp_counter(omp_id+1) + 1
                    list_of_ind_bool(ind_3) = .TRUE.
                end if
            else
                write(*,*) 'some error in 1D list'
            end if
        end do
        !$OMP END PARALLEL DO

        write(*,*) 'test of open_mp_counter'
        do i = 1, threads_max
            write(*,*) i, open_mp_counter(i)
        end do

        ja_val_arr_size = sum(open_mp_counter)
        deallocate( open_mp_counter )
        write(*,*) 'necessary ja and val_arr size = ', ja_val_arr_size

        allocate( list_of_ind_2(ja_val_arr_size, 2) )
        ! allocate( list_of_ind_bool_2(ja_val_arr_size) )
        
        write(*,*) 'Memory study 2: list_of_ind_2: ', ja_val_arr_size
        write(*,*) 'Memory study 2: list_of_ind_2: ', kind(list_of_ind_2(1,1)) * size(list_of_ind_2) , 'bytes', kind(list_of_ind_2(1,1)) * size(list_of_ind_2)/1024.0/1024.0, 'MB'

        ind_temp = 1
        do ind_3 = 1, size_of_list
            if( list_of_ind_bool(ind_3) ) then
                call index_map(Sz_subspace_size, ind_3, ind_Sz_1, ind_Sz_2)
                list_of_ind_2(ind_temp, 1) = ind_Sz_1 !list_of_ind(ind_3, 1)
                list_of_ind_2(ind_temp, 2) = ind_Sz_2 !list_of_ind(ind_3, 2)
                !write(*,*) 'now printed', ind_3, list_of_ind_2(ind_temp, 1), list_of_ind_2(ind_temp, 2)
                ind_temp = ind_temp+1
            end if
        end do

        deallocate(list_of_ind_bool)

        allocate( ia( Sz_subspace_size+1), ja(ja_val_arr_size), val_arr(ja_val_arr_size) )
        !ia and ja inside list_of_ind_2
        
        write(*,*) 'Memory study 3: ia, ja, val_arr: ', Sz_subspace_size+1, ja_val_arr_size, ja_val_arr_size
        write(*,*) 'Memory study 3: ia: ', kind(ia(1)) * size(ia) , 'bytes',  kind(ia(1)) * size(ia)/1024.0/1024.0, 'MB'
        write(*,*) 'Memory study 3: ja: ', kind(ja(1)) * size(ja) , 'bytes',  kind(ja(1)) * size(ja)/1024.0/1024.0, 'MB'
        write(*,*) 'Memory study 3: val_arr: ', kind(val_arr(1)) * size(val_arr) , 'bytes',  kind(val_arr(1)) * size(val_arr)/1024.0/1024.0, 'MB'

        !omp loop
        !asynchronous filling of val_arr!
        !$OMP PARALLEL DO PRIVATE(ind_3, ind_i, ind_j, H_full_ij, omp_id ) SHARED(ja_val_arr_size, hash_Sz, list_of_ind_2, N_spin, J_spin, val_arr, ja) default(none)
        do ind_3 = 1, ja_val_arr_size
            omp_id = omp_get_thread_num()
            ind_i = hash_Sz( list_of_ind_2(ind_3, 1) )
            ind_j = hash_Sz( list_of_ind_2(ind_3, 2) )
            call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
            !write(*,*) 'id ', omp_id , 'ind_3', ind_3, list_of_ind_2(ind_3, 1), list_of_ind_2(ind_3, 2) !  ind_i, ind_j, H_full_ij
            val_arr(ind_3) = real(H_full_ij)
            ja(ind_3) = list_of_ind_2(ind_3, 2)
        end do
        !$OMP END PARALLEL DO

        deallocate(hash_Sz)
        
        ! loop to calculate ia()
        write(*,*) 'Sz_subspace_size', Sz_subspace_size
        ind_temp = 1
        ind_temp_2 = 1
        do ind_3 = 1, Sz_subspace_size-1
            ia(ind_3) = ind_temp
            !write(*,*) 'ind_3, ia', ind_3, ia(ind_3)
            !write(*,*) 'a', list_of_ind_2(ind_temp_2, 1), ind_3
            DO WHILE ( list_of_ind_2(ind_temp_2, 1) == ind_3 )
                ind_temp = ind_temp+1
                ind_temp_2 = ind_temp_2 + 1
                !write(*,*) 'b', ind_temp, ind_temp_2
            END DO
        end do
        ia(Sz_subspace_size) = ind_temp
        !write(*,*) 'ind_3, ia', Sz_subspace_size, ia(Sz_subspace_size)
        ia(Sz_subspace_size +1) = ind_temp+1
        !write(*,*) 'ind_3, ia', Sz_subspace_size +1, ia(Sz_subspace_size +1)
        
        deallocate(list_of_ind_2)
        write(*,*) 'after some deallocations'
        !write(*,*) 'ia: '
        !do ind_3 = 1, Sz_subspace_size+1
        !    write(*,*) ind_3, ia(ind_3)
        !end do

        !write(*,*) 'ja and val_arr:'
        ! do ind_3 = 1, ja_val_arr_size
        !    write(*,*) ind_3, ja(ind_3), val_arr(ind_3)
        !end do

        call calc_timer%stop()
        print *, "----------------------------------------"
        print *, "END CSR filling:"
        call calc_timer%print_elapsed(time_unit%seconds, "seconds")
        call calc_timer%print_elapsed(time_unit%minutes, "minutes")
        call calc_timer%print_elapsed(time_unit%hours, "hours")
        print *, "----------------------------------------"

        call calc_timer%reset()


        call calc_timer%start()
        print *, "----------------------------------------"
        print *, "START FEAST diagonalisation"
        print *, "----------------------------------------"

        !FEAST diagonalization
        write(*,*) 'val_arr size: ', ja_val_arr_size
        write(*,*) 'Sz_subspace_size: ', Sz_subspace_size
        write(*,*) 'size_of_list', size_of_list

        !int 64 has 64 bits
        !double precision has 64 bit

        !n=n     ! Sets the size of the problem
        !a=non_zero_array     ! Array containing the nonzero elements of the upper triangular part of the matrix A
        call feastinit (fpm)  ! function specifying default parameters fpm of FEAST algorithm
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
        
        uplo='U' ! If uplo = 'U', a stores the upper triangular parts of A.
        n = Sz_subspace_size !
        !Intervals for 10 lowest eigenstates for 1D chain of H_XXX with NN hoping
        emin = -0.4465d0 * N_spin + 0.1801d0
        emax = -0.49773d0 * N_spin + 2.10035d0

        write(*,*) "Calculated lower bound: ", emin
        write(*,*) "Calculated upper bound: ", emax

        m0 = 20 !Sz_subspace_size !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n.
        !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax].
        !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        
        allocate( x(n,m0), e(m0), res(m0) )
        write(*,*) 'Memory study 4: n,m0: ', n, m0
        write(*,*) 'Memory study 4: x: ', kind(x(1,1)) * size(x) , 'bytes',  kind(x(1,1)) * size(x)/1024.0/1024.0, 'MB'
        write(*,*) 'Memory study 4: e: ', kind(e(1)) * size(e) , 'bytes',  kind(e(1)) * size(e)/1024.0/1024.0, 'MB'
        write(*,*) 'Memory study 4: res: ', kind(res(1)) * size(res) , 'bytes',  kind(res(1)) * size(res)/1024.0/1024.0, 'MB'
        
        
        write(*,*) 'Before dfeast_scsrev... '
        call dfeast_scsrev(uplo, n, val_arr, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) ' dfeast_scsrev info=', info
        write(*,*) 'After  dfeast_scsrev... '

        if (info /= 0) then
            write(*,*) 'problem with  dfeast_scsrev, info=', info
        end if

        write(*,*) ' dfeast_scsrev eigenvalues found= ', m
        
        ! console print of eigenvals for debug
        do i=1,m
            norm = 0.0
            do j=1, n
                norm = norm + x(j,i)*(x(j,i))
            end do
            write(10,*) i, e(i), norm
        end do

        do i=1,m
            norm = 0.0
            do j=1, n
                norm = norm + x(j,i)*(x(j,i))
            end do
            write(io_unit,*) i, e(i), norm
        end do

        deallocate(val_arr, ia, ja, x, e, res)
        close(10)

        call calc_timer%stop()
        print *, "----------------------------------------"
        print *, "END FEAST diagonalisation:"
        call calc_timer%print_elapsed(time_unit%seconds, "seconds")
        call calc_timer%print_elapsed(time_unit%minutes, "minutes")
        call calc_timer%print_elapsed(time_unit%hours, "hours")
        print *, "----------------------------------------"
        
        end subroutine Hamiltonian_fill_diag_open_mp

end module heisenberg 

program spin_code
    use heisenberg
    use tests_module
    implicit none

    integer :: N_spin, N_spin_max 
    integer (8) :: Sz_subspace_size
    double precision :: J_spin, entropy_value, eigen_value_check, e_up, e_down
    double precision :: start, finish, finish_one_spin, start_one_spin, Sz_choice
    character(len = 12) :: N_spin_char, J_spin_char
    integer, allocatable :: Sz_basis(:), hash_Sz(:), hash(:), indices_Sz_basis_sorted(:), basis_vector(:,:), basis_rho_target(:,:), new_basis(:,:)
    double precision, allocatable :: target_sz(:), eigen_values(:), eigen_vectors(:,:), rho_reduced(:,:)

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
    call omp_mkl_small_test()
    ! call test_permutation_H_for_4_sites()
    ! call sparse_zfeast_test()
    call sparse_dfeast_test()

    write(*,*) '------- END test modules -------'

    !N_spin = 4
    !J_spin = 1 

    call cpu_time(start)
    write(*,*) '------- START Heisenberg Program -------'
    write(*,*) 'Calculation of Heisenberg chain for N =', N_spin , 'and J = ', J_spin 
    write(*,*) ' '

    N_spin_max = 2**N_spin 
    ! TO DO: 
    ! adjust the choose of the sz choice based one the N_spin automatically! 
    ! write the output to the log file 
    ! create result file
    !Sz_choice = 0.0d0 !integer counted from Sz_max (in a sense that Sz_choice = 1 means eg for N_spin=4, Sz_max=2 Sz=2)
    !Choose always biggest subspace of the Hamiltonian
    if (N_spin%2 == 0.0d0) then
        Sz_choice = 0.0d0 
    else 
        Sz_choice = 0.5d0  
    endif 
    
    call Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)

    call Hamiltonian_fill_diag_open_mp(N_spin, J_spin, Sz_subspace_size, hash_Sz)

    write(*,*) " "
    write(*,*) "Program executed with success"
    call cpu_time(finish)
    print '("Time = ",f," seconds.")',finish-start
    write(*,*) '------- END Heisenberg Program -------'

    

end program spin_code