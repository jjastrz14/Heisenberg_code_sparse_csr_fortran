include "mkl_vsl.f90" !necessary for intel mkl random numbers

module spin_systems
    contains
    subroutine binomialCoefficient(n, k, C)
        integer, intent(in)  :: n, k
        integer(8), intent(out) :: C
        !use, intrinsic :: iso_fortran_env, only: real64, int64
        !integer(int64) :: x
        integer :: i
        integer :: C_temp

        C_temp = 1

        if (k < 0 .or. k > n) then
        C_temp = 0
        else
        do i = 1, min(k, n-k)
            C_temp = C_temp * (n - i + 1) / i
            if (C_temp.GT.2147483647) then
                write(*,*) 'possible integer 32bit 2^32 = 2147483647 overflow'
                error stop
            end if
        end do
        end if

        C = C_temp
    end subroutine binomialCoefficient

    subroutine index_map(Sz_subspace_size, ind_1D_target, ind_i, ind_j) 
        implicit none 
        integer(8), intent(in) :: Sz_subspace_size, ind_1D_target
        integer, intent(out) :: ind_i, ind_j
        integer(8) :: ind_3
        integer(8) :: ind_Sz_1, ind_Sz_2, ind_tnp, loop_occ, N

        inD_tmp = ind_1D_target
        N = Sz_subspace_size
        loop_occ = 0 

        do while (int_tmp .GT. 0)
            ind_tmp = ind_tmp - N
            N = N - 1
            loop_occ = loop_occ + 1
        end do 
        ind_i = loop_occ
        ind_j = ind_i + ind_tmp + N 

    end subroutine index_map


    subroutine Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)
        implicit none
        integer, intent(in) :: N_spin
        double precision, intent(in) :: Sz_choice
        integer(8), intent(out) :: Sz_subspace_size !beacuse of overflow for N_spin >= 18
        integer, allocatable, intent(out) :: hash_Sz(:) !try to allocate it as small as possible to save memory
        
        integer :: N_spin_max, i, j, ind_2, Sz_choice_ind
        double precision :: Sz
        logical :: bool
        
        !double precision, allocatable   :: Sz_basis(:)
        !info
        !
        !CHARACTER*1 :: id

        N_spin_max = 2**N_spin
        
        !from all N_spin_max we should be able to predict degeneracy of restricted S_z
        !1) we need to know Sz_choice_ind from Sz_choice and N_spin_max
        !Sz_choice = N_spin_max/2      -> Sz_choice_ind = 0
        !Sz_choice = N_spin_max/2-1    -> Sz_choice_ind = 1
        !Sz_choice = N_spin_max/2-2    -> Sz_choice_ind = 2
        write(*,*) 'choosen spin= ', Sz_choice
        
        Sz_choice_ind = int( N_spin/2.0d0 - Sz_choice )
        write(*,*) 'calculated integer Sz_choice_ind = ', Sz_choice_ind
        
        !2) Sz_subspace_size = 4 here Newton symbol to be implemented = ( N_spin_max over Sz_choice_ind )
        
        !testing binomialCoefficient(n, k, C)
        call binomialCoefficient(N_spin, Sz_choice_ind, Sz_subspace_size)
        
        !Sz_subspace_size = 4 ! standard factorial calcualated recursively/ logartihmically
        write(*,*) 'Sz_subspace_size = ', Sz_subspace_size
        
        !Sz_subspace_size = 1
        
        allocate ( hash_Sz(Sz_subspace_size) )

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
        
        !open (unit=40, file="hash_Sz_basis.dat", recl=512)
        !write(*,*) 'Summary of hash_Sz basis: '
        !do i=1,ind_2-1
        !    write(* ,*) i, hash_Sz(i)
        !    write(40,*) i, hash_Sz(i)
        !end do
        !close(40)
    end subroutine Sz_subspace_choice


    subroutine H_XXX_subspace_fill_and_diag_fast_time_window(N_spin, J_spin, Sz_choice, Sz_subspace_size, hash_Sz, e, x, m)
        use omp_lib
        use mkl_vsl
        implicit none
        !character(len = 12), intent(in) :: N_spin_char
        integer, dimension(8) :: start_csr, finish_csr, start_diag, finish_diag
        integer, intent(in) :: N_spin
        integer (8), intent(in) :: Sz_subspace_size
        double precision, intent(in) :: Sz_choice
        integer, allocatable :: hash_Sz(:), list_of_ind_2(:,:), ia(:), ja(:), open_mp_counter(:) !list_of_ind(:,:) changed to 2 bytes below
        logical, allocatable :: list_of_ind_bool(:)
        double precision, intent(in) :: J_spin
        integer :: i, j, ind_i, ind_j, ind_ja, ind_temp, ind_temp_2, N_spin_max, size_of_1D_list_for_sweep, &
                    ja_val_arr_size, omp_id, threads_max
        integer(4), allocatable :: list_of_ind(:,:) ! Signed integer value from -32,768 to 32,767
        integer(8) :: ind_Sz_1, ind_Sz_2, ind_3
        integer(8) :: size_of_list !8 byte = 64 bit Signed integer value from -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
        double precision, allocatable :: H_full(:,:), val_arr(:)
        double precision :: norm, H_full_ij, element_value_cutoff, e_min_reg, e_max_reg
        double precision, allocatable, intent(out) :: x(:,:)
        integer :: fpm(128)
        double precision :: emin, emax
        integer :: m0
        double precision :: epsout, t1, t2
        integer :: loop
        double precision, allocatable, intent(out) :: e(:)
        integer, intent(out) :: m
        double precision, allocatable :: res(:)
        
        CHARACTER*1 :: jobz, range, uplo
        
        integer :: il, iu, ldz, liwork, lwork, info, m_eig, n

        ! lets rewrite our matlab code
        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        call date_and_time(VALUES=start_csr)
        t1 = omp_get_wtime()

        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size - large integer possibility'
        write(*,*) 'test for N_spin_max', N_spin_max
        
        ! 1D list for fast sweep
        size_of_list = Sz_subspace_size*(Sz_subspace_size + 1  ) / 2

        write(*,*) 'Before the lists allocation for 1d sweep'
        write(*,*) 'size_of_list = ', size_of_list
        allocate( list_of_ind(size_of_list, 2) )
        write(*,*) 'After the list_of_ind allocation for 1d sweep'
        allocate( list_of_ind_bool(size_of_list) )
        write(*,*) 'After the list_of_ind bool allocation for 1d sweep'

        list_of_ind_bool = .FALSE.
        ind_3 = 1
        do ind_Sz_1 = 1        , Sz_subspace_size
        do ind_Sz_2 = ind_Sz_1 , Sz_subspace_size
            list_of_ind(ind_3,1) = ind_Sz_1
            list_of_ind(ind_3,2) = ind_Sz_2
            ind_3 = ind_3 + 1
        end do
        end do
        
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
        !$OMP PARALLEL DO PRIVATE(ind_3, omp_id, ind_Sz_1, ind_Sz_2, ind_i, ind_j, H_full_ij) &
                        SHARED(size_of_list, list_of_ind, open_mp_counter, list_of_ind_bool, hash_Sz, N_spin, J_spin, element_value_cutoff) &
                        default(none)
                        
        do ind_3 = 1, size_of_list
            omp_id = omp_get_thread_num()
            ind_Sz_1 = list_of_ind(ind_3,1)
            ind_Sz_2 = list_of_ind(ind_3,2)
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
            !if (list_of_ind_bool(ind_3)) then
            !    open_mp_counter(omp_id+1) = open_mp_counter(omp_id+1) + 1
            !end if
        end do
        !$OMP END PARALLEL DO 
        
        write(*,*) 'test of open_mp_counter'
        do i = 1, threads_max
            write(*,*) i, open_mp_counter(i)
        end do
        
        ja_val_arr_size = sum(open_mp_counter)
        write(*,*) 'necessary ja and val_arr size = ', ja_val_arr_size
        
        
        !do ind_3 = 1, size_of_list
            !write(*,*) ind_3, list_of_ind(ind_3,1), list_of_ind(ind_3,2), list_of_ind_bool(ind_3)                    
        !end do
        write(*,*) 'Before the lists allocation for 2nd sweep'
        allocate( list_of_ind_2(ja_val_arr_size, 2) )
    ! allocate( list_of_ind_bool_2(ja_val_arr_size) )
        
        ind_temp = 1
        do ind_3 = 1, size_of_list
            if( list_of_ind_bool(ind_3) ) then
                list_of_ind_2(ind_temp, 1) = list_of_ind(ind_3, 1)
                list_of_ind_2(ind_temp, 2) = list_of_ind(ind_3, 2)
                !write(*,*) ind_3, list_of_ind_2(ind_temp, 1), list_of_ind_2(ind_temp, 2)
                ind_temp = ind_temp+1
            end if
        end do
        
        deallocate (list_of_ind, list_of_ind_bool)
        write(*,*) 'Before the csr vectors allocation'
        allocate( ia( Sz_subspace_size+1), ja(ja_val_arr_size), val_arr(ja_val_arr_size) )
        !ia and ja inside list_of_ind_2
        
        !asynchronous filling of val_arr!
        !$OMP PARALLEL DO PRIVATE(ind_ja, ind_i, ind_j, H_full_ij, omp_id) &
                        SHARED(ja_val_arr_size, hash_Sz, list_of_ind_2, N_spin, J_spin, val_arr, ja) &
                        default(none)
        do ind_ja = 1, ja_val_arr_size      
            omp_id = omp_get_thread_num()
            ind_i = hash_Sz( list_of_ind_2(ind_ja, 1) )
            ind_j = hash_Sz( list_of_ind_2(ind_ja, 2) )
            call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
            !write(*,*) 'id ', omp_id , 'ind_3', ind_3, ind_i, ind_j, H_full_ij
            val_arr(ind_ja) = H_full_ij
            ja(ind_ja) = list_of_ind_2(ind_ja, 2)
        end do
        !$OMP END PARALLEL DO 
        
        ! loop to calculate ia()
        ind_temp = 1
        ind_temp_2 = 1
        do ind_3 = 1, Sz_subspace_size                
            ia(ind_3) = ind_temp
            DO WHILE ( list_of_ind_2(ind_temp_2, 1) == ind_3 ) 
                ind_temp = ind_temp+1
                ind_temp_2 = ind_temp_2 + 1
            END DO 
        end do
        ia(Sz_subspace_size +1) = ind_temp
        
        !write(*,*) 'ia: '
        !do ind_3 = 1, Sz_subspace_size+1
        !    write(*,*) ind_3, ia(ind_3)    
        !end do
        
        !write(*,*) 'ja and val_arr:'
        !do ind_3 = 1, ja_val_arr_size
        !    write(*,*) ind_3, ja(ind_3), val_arr(ind_3)
        !end do
        call date_and_time(VALUES=finish_csr) !end time counter
        t2 = omp_get_wtime()

        !old example for date_and_time
        !write(*,*) 'Csr filling time = ', &
            !(((finish_csr(5) - start_csr(5))*3600 + &
            !(finish_csr(6) - start_csr(6))*60 + &
            !(finish_csr(7) - start_csr(7)))*1000 + &
            !(finish_csr(8) - start_csr(8)))/1000, "seconds."

        print '("Csr filling time = ",f," seconds.")',t2-t1
        write(3,*) "Csr filling time = ", t2-t1 ," seconds."

        !Calculation of energy bounds for the chain problem according to the previous estimations: 
        e_min_reg = -0.4465d0 * N_spin + 0.1801d0
        e_max_reg = -0.49773d0 * N_spin + 2.10035d0

        write(*,*) "Calculated lower bound: ", e_min_reg
        write(*,*) "Calculated upper bound: ", e_max_reg
        !FEAST diagonalization

        !n=n     ! Sets the size of the problem
        !a=non_zero_array     ! Array containing the nonzero elements of the upper triangular part of the matrix A
        call feastinit (fpm)  ! function specifying default parameters fpm of FEAST algorithm
        fpm(1) = 1
        fpm(2) = 12 !can be more, can be less
        fpm(3) = 8 !eps 10^-fpm(3)
        fpm(4) = 20 ! max number of feast loops
        fpm(5) = 0 !initial subspace
        fpm(6) = 0! stopping criterion
        fpm(7) = 5 !Error trace sigle prec stop crit
        fpm(14) = 0 ! standard use of feast
        fpm(27) = 1 !check input matrices
        fpm(28) = 1 !check if B is positive definite?         
        uplo='U' ! If uplo = 'U', a stores the upper triangular parts of A.
        emin = e_min_reg! The lower ... &
        emax = e_max_reg  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = 12 !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0?n. 
        !Set m0 ? m where m is the total number of eigenvalues located in the interval [emin, emax]. 
        !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = Sz_subspace_size
        allocate( x(n,m0), e(m0), res(m0) )
        !write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
        write(*,*) 'Before dfeast_scsrev... '

        call date_and_time(VALUES=start_diag)
        t1 = omp_get_wtime()
        call dfeast_scsrev(uplo, n, val_arr, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
        t2 = omp_get_wtime()
        call date_and_time(VALUES=finish_diag)
        
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) ' dfeast_scsrev info=', info
        write(*,*) 'After  dfeast_scsrev... '
        
        if (info /= 0) then
        write(*,*) 'problem with  dfeast_scsrev, info=', info
        end if 

        print '("Diagonalisation time = ",f," seconds.")',t2-t1
        write(3, *) "Diagonalisation time = ", t2-t1 ," seconds."

        !file_name1 = 'Eigenvalues_results_' // trim(adjustl(N_spin_char)) // '_feast.dat'
        !open (unit=1, file= trim(file_name1), recl=512)
        !write(1,*), "#Eigenvalue     Spin_z   norm" 

        !write(*,*) ' dfeast_scsrev eigenvalues found= ', m
        !do i = 1, m
         !   write(*,*) i, e(i)
        !end do

        !write(*,*) ' dfeast_scsrev eigenvec norm:'
        !do i=1,m
        !    norm = 0.0d0
         !   do j=1, n
         !       norm = norm + x(j,i)*(x(j,i))
         !   end do
         !   write(*,*) i, norm
          !  write(1,*) e(i), ',' , Sz_choice, ',', norm
            
        !end do

        !close(1)      
        !deallocate( val_arr, ia, ja, x, e, res)
        deallocate( val_arr, ia, ja, res)
        deallocate( list_of_ind_2 )
        deallocate( open_mp_counter )

    end subroutine H_XXX_subspace_fill_and_diag_fast_time_window


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

end module spin_systems

program spin_code
    use spin_systems
    use omp_lib
    implicit none

    integer :: N_spin, N_spin_test, N_spin_max, no_of_nonzero, size_of_sub_A, size_of_sub_B, i,j, k, target_sz_spin, number_of_eigen, size_of_eigenvalues, m
    integer (8) :: Sz_subspace_size
    double precision :: J_spin, entropy_value, eigen_value_check, e_up, e_down, norm
    double precision :: start, finish, finish_one_spin, start_one_spin, start_csr, finish_csr, start_diag, finish_diag, Sz_choice
    character(len = 12) :: N_spin_char, J_spin_char
    character(len=53) :: dec2bin, spin_basis, file_name1, file_name2, file_name3
    integer, allocatable :: Sz_basis(:), hash_Sz(:), hash(:), indices_Sz_basis_sorted(:), basis_vector(:,:), basis_rho_target(:,:), new_basis(:,:)
    double precision, allocatable :: target_sz(:), eigen_values(:), eigen_vectors(:,:), rho_reduced(:,:), e(:), x(:,:)

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

    !11.03/2024 - powrót do kodu :)

    If(command_argument_count().NE.2) Then
        write(*,*)'Error, Only N (integer) and J (double precision) is required, program stopped two many or parameter is missing'
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
    ! call sparse_zfeast_test()
    call sparse_dfeast_test()

    start = omp_get_wtime()
    write(*,*) '------- START Heisenberg Program -------'
    N_spin_max = 2**N_spin 
    allocate(target_sz(N_spin-1))

    !file_name1 = 'Test_Eigenvalues_results_feast.dat'
    file_name2 = 'Eigenvalues_results_' // trim(adjustl(N_spin_char)) // '_feast.dat'
    file_name3 = 'Cpu_Time_details_' // trim(adjustl(N_spin_char)) // '.dat'

    !open (unit=1, file= trim(file_name1), recl=512)
    open (unit=2, file= trim(file_name2), recl=512)
    open (unit=3, file= trim(file_name3), recl=512)
    !write(1,*), "#Eigenvalue     Spin_z   norm" 
    write(2,*), "#Eigenvalue     Spin_z   norm" 
    write(3,*), "#Time needed for each process"
    
    ! generate a target_sz array storing all possible Sz values for N_spin 
    
    do i = 1, N_spin - 1
            target_sz(i) = N_spin/2.0d0 - i
    end do 


    print *, "All combination of S_z spin", target_sz
    write(*,*) " "
    
    !write(*,*) '------- START Heisenberg diagonalisation test -------'
    !N_spin_test = 4
    !Sz_choice = 0.0d0 ! integer counted from Sz_max (in a sense that Sz_choice = 1 means eg for N_spin=4, Sz_max=2 Sz=2)
    !call Sz_subspace_choice(N_spin_test, Sz_choice, hash_Sz, Sz_subspace_size)
    ! above we confirmed Sz_choice works
    !write(3,*), "#Heisenberg diagonalisation test: N = 4"
    !call H_XXX_subspace_fill_and_diag_fast(N_spin, J_spin, Sz_choice, Sz_subspace_size, hash_Sz, e, x, m) !here only csr3

    !write(*,*) ' dfeast_scsrev eigenvalues found= ', m
    !do i = 1, m
    !    write(*,*) i, e(i)
    !end do

    !write(*,*) ' dfeast_scsrev eigenvec norm:'
    !do i=1,m
    !    norm = 0.0d0
    !    do j=1, Sz_subspace_size
    !        norm = norm + x(j,i)*(x(j,i))
    !    end do
    !    write(*,*) i, norm
    !   write(1,*) e(i), ',' , Sz_choice, ',', norm
        
    !end do

    !write(*,*) '------- END Heisenberg diagonalisation test -------'

    write(*,*) '------- Start Heisenberg diagonalisation loop -------'

    !here write if describing the situation if the size of the matrix is an integer ;)

    !do k = 1, size(target_sz)
    !Sz_choice = target_sz(k)
    if (mod(N_spin, 2) == 0) then
        ! This is the case for even numbers
        Sz_choice = 0.0d0
    else
        ! This is the case for odd numbers
        Sz_choice = 0.5d0
    endif

    write(*,*) 'Start loop for Sz =', Sz_choice
    write(3,*), "#Heisenberg diagonalisation: N = " , N_spin_char

    call Sz_subspace_choice(N_spin, Sz_choice, hash_Sz, Sz_subspace_size)
    call H_XXX_subspace_fill_and_diag_fast_time_window(N_spin, J_spin, Sz_choice, Sz_subspace_size, hash_Sz, e, x, m)

    write(*,*) ' dfeast_scsrev eigenvalues found= ', m
    write(*,*) ' dfeast_scsrev eigenvec norm:'
    do i=1,m
        norm = 0.0d0
        do j=1, Sz_subspace_size
            norm = norm + x(j,i)*(x(j,i))
        end do
        write(2,*) e(i), ',' , Sz_choice, ',', norm
    end do 

    deallocate(e, x, hash_Sz)

    !end do

    deallocate(target_sz)
    write(*,*) '------- END Heisenberg diagonalisation loop -------'  
    write(*,*) " "
    write(*,*) "Program executed with success"
    finish = omp_get_wtime()
    print '("Final time = ",f," seconds.")',finish-start
    write(3, *) "Final time = ", finish-start ," seconds." 
    close(1) 
    close(2)  
    close(3)
    write(*,*) '------- END Heisenberg Program -------'


end program spin_code