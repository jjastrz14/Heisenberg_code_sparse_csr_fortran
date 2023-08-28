module spin_systems
    contains

    subroutine H_create_basis_sz(N_spin, indices_Sz_basis_sorted)
        !> Returns a binary representation of basis(i) as a string with at least N_spin bits. 
        !> abs(d) < 2^52
        ! dec2bin(d, n)
        use math_functions
        implicit none

        integer, intent(in) :: N_spin
        double precision :: Sz
        character(len=53), allocatable :: dec2bin(:)
        double precision, allocatable   :: basis_sz(:), Sz_basis(:), Sz_basis_target(:)
        integer :: N_spin_max, i, j, info, N, target_sz
        integer, allocatable, intent(out) :: indices_Sz_basis_sorted(:)
        integer, allocatable :: basis(:), hash(:)
        logical :: bool
        CHARACTER*1 :: id
        character(len=53)             :: tmp
        integer                       :: n_
        character(len=8)              :: f

        N_spin_max = 2**N_spin

        allocate (Sz_basis(N_spin_max), indices_Sz_basis_sorted(N_spin_max))
        allocate (basis(N_spin_max))
        allocate (dec2bin(N_spin_max))
        !allocate (index_array(N_spin_max))

        ! only for visualization of the basis by 0 and 1 combinations
        if (N_spin <= 10) then
            do i = 1 , N_spin_max
                basis(i) = i
            end do

            !now we generate the basis +1/2 -> 0, -1/2 -> 1
            do i = 1, N_spin_max
            n_ = min(N_spin, 53)
            write(f,'(i2)') n_
                f = '(B' // trim(adjustl(f)) // '.' // trim(adjustl(f)) // ')'
                write(tmp,f) basis(i) - 1
                dec2bin(i) = trim(adjustl(tmp))
            end do

            write(*,*) 'Binary basis: ' 
            write(*,*) dec2bin
        end if

        ! Maciek version of basis generation with btest
        ! we can use e.g. btest function, which returns info if given bit is 0 (false) or 1 (true) 
        !bool_T/F=btest(number, bit_number)
        ! 0 == false => Sz=Sz+1/2, 1= true> = Sz=Sz-1/2
        do i = 0 , N_spin_max-1
            !basis(i) = i !here we want Sz from this i
            ! write(*,*) 'basis vector i =', i+1
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
            Sz_basis(i+1) = Sz
            ! write(*,*) 'for basis vector i =', i, 'Sz = ', Sz_basis(i+1)
        end do

        write(*,*) 'Summary of Sz basis: '
        do i=1,N_spin_max
            write(*,*) Sz_basis(i)
        end do

        !sorting
        id = 'D'
        !dlasrt
        call dlapst(id, N_spin_max, Sz_basis, indices_Sz_basis_sorted, info ) !scalapack quicksort

        !subroutine sorts basis_sz array and creates array for permutation matrix
        ! call indexArrayReal(N_spin_max, basis_sz, index_array)

        write(*,*) 'Indicies of Sorted Sz basis: '
        write(*,*) indices_Sz_basis_sorted

        write(*,*) 'Sorted Sz basis: '
        write(*,*) Sz_basis(indices_Sz_basis_sorted)

        deallocate(Sz_basis)

    end subroutine H_create_basis_sz


    subroutine H_create_basis_sz_with_target(N_spin, hash, target_sz)
        !> Returns a binary representation of basis(i) as a string with at least N_spin bits. 
        !> abs(d) < 2^52
        ! dec2bin(d, n)
        use math_functions
        implicit none

        integer, intent(in) :: N_spin, target_sz
        double precision :: Sz
        character(len=53), allocatable :: dec2bin(:)
        double precision, allocatable   :: basis_sz(:), Sz_basis(:), S_z_target(:)
        integer :: N_spin_max, i, j, info, N
        integer, allocatable, intent(out) ::  hash(:)
        integer, allocatable :: basis(:), indices_Sz_basis_sorted(:)
        logical :: bool
        CHARACTER*1 :: id
        character(len=53)             :: tmp
        integer                       :: n_
        character(len=8)              :: f

        N_spin_max = 2**N_spin

        allocate (Sz_basis(N_spin_max))
        allocate (basis(N_spin_max))
        allocate (dec2bin(N_spin_max))
        allocate (indices_Sz_basis_sorted(N_spin_max))

        ! only for visualization of the basis by 0 and 1 combinations
        if (N_spin <= 10) then
            do i = 1 , N_spin_max
                basis(i) = i
            end do

            !now we generate the basis +1/2 -> 0, -1/2 -> 1
            do i = 1, N_spin_max
            n_ = min(N_spin, 53)
            write(f,'(i2)') n_
                f = '(B' // trim(adjustl(f)) // '.' // trim(adjustl(f)) // ')'
                write(tmp,f) basis(i) - 1
                dec2bin(i) = trim(adjustl(tmp))
            end do

            write(*,*) 'Binary basis: ' 
            write(*,*) dec2bin
        end if

        ! Maciek version of basis generation with btest
        ! we can use e.g. btest function, which returns info if given bit is 0 (false) or 1 (true) 
        !bool_T/F=btest(number, bit_number)
        ! 0 == false => Sz=Sz+1/2, 1= true> = Sz=Sz-1/2
        do i = 0 , N_spin_max-1
            !basis(i) = i !here we want Sz from this i
            ! write(*,*) 'basis vector i =', i+1
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
            Sz_basis(i+1) = Sz
            ! write(*,*) 'for basis vector i =', i, 'Sz = ', Sz_basis(i+1)
        end do

        write(*,*) 'Summary of Sz basis: '
        do i=1,N_spin_max
            write(*,*) Sz_basis(i)
        end do

        !sorting
        id = 'D'
        !dlasrt
        call dlapst(id, N_spin_max, Sz_basis, indices_Sz_basis_sorted, info ) !scalapack quicksort

        !subroutine sorts basis_sz array and creates array for permutation matrix
        ! call indexArrayReal(N_spin_max, basis_sz, index_array)

        write(*,*) 'Sorted Sz basis: '
        write(*,*) Sz_basis(indices_Sz_basis_sorted)

        allocate(hash(N_spin_max), S_z_target(N_spin_max))

        S_z_target = Sz_basis(indices_Sz_basis_sorted)

        write(*,*) 'Sorted Sz target basis: '
        write(*,*) S_z_target

        N = 1
        do i = 1, N_spin_max
            if (S_z_target(i) == target_sz) then
                !Sz_basis_target(i) = i
                !Sz_basis_target(N) = i
                hash(i) = N  
                N = N + 1
            else 
                hash(i) = -1 
            endif 
        end do 

        write(*,*) 'Hash:'
        write(*,*) hash

        deallocate(Sz_basis, indices_Sz_basis_sorted, S_z_target)

    end subroutine H_create_basis_sz_with_target

    subroutine Create_permutation_matrix(N_spin_max, array_index, p_matrix)
        implicit none

        integer, intent(in) :: N_spin_max
        integer, intent(out) :: p_matrix(N_spin_max, N_spin_max)
        integer, intent(in)  :: array_index(N_spin_max)
        integer :: i, j

        do i = 1, N_spin_max
            p_matrix(i,array_index(i)) = 1
        end do
        
    end subroutine Create_permutation_matrix


    subroutine Create_permutation_matrix_csr(N_spin_max, array_index, no_of_nonzero)
        implicit none

        integer, intent(in) :: N_spin_max
        integer, intent(in)  :: array_index(N_spin_max)
        integer :: i, j, ind_i, ind_j, counter, ind_ia
        double precision :: P_full_ij 
        integer, intent(out):: no_of_nonzero

        ! function for creating a permutation matrix in csr format 

        open (unit=40, file="ia_per.dat", recl=512)
        open (unit=41, file="ja_per.dat", recl=512)
        open (unit=42, file="values_array_per.dat", recl=512)

        do i = 1, N_spin_max

            write(40, *) i
            write(41, *) array_index(i)
            write(42, *) 1.0

        end do
        write(40, *) i
    
        no_of_nonzero = N_spin_max
        write(*,*) 'number of non-zero elements = ', no_of_nonzero
        
        close(40)
        close(41)
        close(42)

    end subroutine Create_permutation_matrix_csr

    subroutine CSR_matrix_multiplication_for_3_matrices(N_spin, J_spin, index_array)
        use mkl_spblas
        use ISO_C_BINDING
        implicit none 
                    
        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer :: N_spin_max, i, j, ja_temp, ia_temp, no_of_nonzero_p_matrix, no_of_nonzero_h_matrix, stat_permutation, stat_H_csr, stat_csr_mutliplication, info
        double precision :: val_arr_temp

        double precision, allocatable :: values_array(:), values_array_per(:)
        integer, allocatable :: ia(:), ja(:), ia_per(:), ja_per(:)
        integer, allocatable, intent(in) :: index_array(:)

        integer :: nCol,nrowsD, ncolsD
        
        !export from internal sparse to CSR
        integer(C_INT) :: indexing
        type(C_PTR)    :: rowsD_start, rowsD_end, colD_indx, Dvalues
        integer   , POINTER :: rowsD_start_f(:), rowsD_end_f(:), colD_indx_f(:)
        double precision, POINTER :: Dvalues_f(:)

        !BLAS types
        type(sparse_matrix_t) :: p_matrix_csr
        type(sparse_matrix_t) :: H_matrix_csr
        type(sparse_matrix_t) :: H_matrix_permuted
        type(matrix_descr) :: descrB, descrP

        N_spin_max = 2**N_spin
        !GENERAL CSR3 FORMAT
        call Create_permutation_matrix_csr(N_spin_max, index_array, no_of_nonzero_p_matrix)

        allocate(values_array_per(no_of_nonzero_p_matrix), ia_per(N_spin_max+1), ja_per(no_of_nonzero_p_matrix))
        
        !reading of exsisting files
        open (unit=40, file="ia_per.dat", recl=512)
        open (unit=41, file="ja_per.dat", recl=512)
        open (unit=42, file="values_array_per.dat", recl=512)
        
        do i=1, N_spin_max+1           
            read(40, *) ia_temp
            ia_per(i) = ia_temp
        end do
            
        do i=1, no_of_nonzero_p_matrix
            read(42,*) val_arr_temp
            values_array_per(i) = val_arr_temp
            read(41, *) ja_temp
            ja_per(i) = ja_temp
        end do
        
        close(40)
        close(41)
        close(42)

        !stat = mkl_sparse_d_create_csr(A, indexing, rows, cols, rows_start, rows_end, col_indx, values)

        !rewriting permutation matrix to intel mkl internal csr format !BLAS POINTERS TO CSR in MKL/BLAS
        stat_permutation = mkl_sparse_d_create_csr(p_matrix_csr, SPARSE_INDEX_BASE_ONE, N_spin_max, N_spin_max , ia_per(1:N_spin_max), ia_per(2:N_spin_max+1), ja_per, values_array_per)
        print *, "stat CSR Permutation matrix create = ", stat_permutation

        !   Create matrix descriptor
        descrP % TYPE = SPARSE_MATRIX_TYPE_GENERAL

        !   Analyze sparse matrix; chose proper kernels and workload balancing strategy
        info = MKL_SPARSE_OPTIMIZE(p_matrix_csr)

        print *, "optimization of CSR Permutation matrix = ", info

        call H_XXX_feast_vec_fill(N_spin, J_spin, no_of_nonzero_h_matrix)

        allocate( values_array(no_of_nonzero_h_matrix), ia(N_spin_max+1), ja(no_of_nonzero_h_matrix) )
        
        !reading of exsisting files
        open (unit=60, file="ia_h.dat", recl=512)
        open (unit=61, file="ja_h.dat", recl=512)
        open (unit=62, file="values_array_h.dat", recl=512)
        
        do i=1, N_spin_max+1           
            read(60, *) ia_temp
            ia(i) = ia_temp
        end do
            
        do i=1, no_of_nonzero_h_matrix
            read(62,*) val_arr_temp
            values_array(i) = val_arr_temp
            read(61, *) ja_temp
            ja(i) = ja_temp
        end do
        
        close(60)
        close(61)
        close(62)

        !rewriting H matrix to intel mkl internal csr format
        stat_H_csr = mkl_sparse_d_create_csr(H_matrix_csr, SPARSE_INDEX_BASE_ONE, N_spin_max, N_spin_max , ia(1:N_spin_max), ia(2:N_spin_max+1), ja, values_array)
        print *, "stat CSR H matrix create = ", stat_H_csr

        descrB%type = SPARSE_MATRIX_TYPE_SYMMETRIC
        descrB%mode = SPARSE_FILL_MODE_UPPER
        descrB % DIAG = SPARSE_DIAG_NON_UNIT !this shouldn't be used but it is 

        ! Analyze sparse matrix; chose proper kernels and workload balancing strategy
        info = MKL_SPARSE_OPTIMIZE(H_matrix_csr)

        print *, "optimization of H matrix = ", info

        !Doing a product of three csr matrices P @ H @ P.T
        stat_csr_mutliplication = mkl_sparse_sypr(SPARSE_OPERATION_NON_TRANSPOSE, p_matrix_csr, H_matrix_csr, descrB, H_matrix_permuted, SPARSE_STAGE_FULL_MULT)

        print *, "stat P @ H matrix @ P.T  = ", stat_csr_mutliplication

        ! export the Hamiltonian in CSR format matrix from the internal representation
        info = mkl_sparse_d_export_csr(H_matrix_permuted, indexing, nrowsD, ncolsD, rowsD_start, rowsD_end, colD_indx, Dvalues)
        
        print *, "Export H permuted to csr3  = ", info
        info = mkl_sparse_order(H_matrix_permuted)
        print *, "Ordering H permuted  = ", info
    
    !   Converting C into Fortran pointers   to czy na pewno potrzebne?
        call C_F_POINTER(rowsD_start, rowsD_start_f, [nrowsD])
        call C_F_POINTER(rowsD_end  , rowsD_end_f  , [nrowsD]) 
        call C_F_POINTER(colD_indx  , colD_indx_f  , [rowsD_end_f(nrowsD)-indexing])
        call C_F_POINTER(Dvalues    , Dvalues_f    , [rowsD_end_f(nrowsD)-indexing])


        write(*,*) 'values array: '
        write(*,*) Dvalues_f

        write(*,*) 'ia: '
        write(*,*) rowsD_start_f ! powinienes polaczyc rowsD_start_f + rowsD_end_f

        write(*,*) 'ja: '
        write(*,*)  colD_indx_f


        deallocate(values_array_per, ia_per, ja_per)
        deallocate(values_array, ia, ja)

        !   Release internal representation of CSR matrix
        info = MKL_SPARSE_DESTROY(p_matrix_csr)
        info = MKL_SPARSE_DESTROY(H_matrix_csr)
        info = MKL_SPARSE_DESTROY(H_matrix_permuted)

    end subroutine CSR_matrix_multiplication_for_3_matrices

    subroutine H_XXX_full()
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none
        
        integer :: i, j, N_spin, J_spin, N_spin_max
        double complex, allocatable :: H_full(:,:), Sp_site_ind_Sm_site_ind_p1(:,:), Sm_site_ind_Sp_site_ind_p1(:,:), Sz_site_ind_Sz_site_ind_p1(:,:)  
        integer :: site_ind, m, n, p, q, r, s
        integer :: ind_i, ind_j, a_i, a_j, b_i, b_j, c_i, c_j
        double precision :: A_val, B_val, C_val, result_ij, norm
 
        CHARACTER*1 :: jobz, range, uplo
        integer :: il, iu, ldz, liwork, lwork, lrwork, info, m_eig
        double precision :: vl, vu, abstol
        double complex, allocatable :: work_complex(:), eigen_vec(:,:)
        integer, allocatable :: iwork(:), isuppz(:)
        double precision, allocatable :: rwork(:), w(:)
        
        open (unit=10, file="H_check.dat", recl=512)
        ! lets rewrite our matlab code
        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        
        call omp_mkl_small_test()

        N_spin = 4
        J_spin = 1
        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
    
    
        allocate( H_full(N_spin_max, N_spin_max), Sp_site_ind_Sm_site_ind_p1(N_spin_max, N_spin_max), Sm_site_ind_Sp_site_ind_p1(N_spin_max, N_spin_max), Sz_site_ind_Sz_site_ind_p1(N_spin_max, N_spin_max) )
    
        H_full = dcmplx(0.0d0,0.0d0)
       
        write(*,*) 'H_full(1,1)', H_full(1,1)

    
        do site_ind = 1, (N_spin-1)
            m = 2**(site_ind-1)
            n = m
            p = 4 ! N_spin independent because this referes to S+S- operator for spin 1/2 
            q = 4
            r = 2**(N_spin-(site_ind+1))
            s = r
            !write(*,*) site_ind, m, n, p, q, r, s
        
            Sp_site_ind_Sm_site_ind_p1 = dcmplx(0.0d0,0.0d0)
            Sm_site_ind_Sp_site_ind_p1 = dcmplx(0.0d0,0.0d0)
            Sz_site_ind_Sz_site_ind_p1 = dcmplx(0.0d0,0.0d0)
            
            do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max
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
                result_ij = A_val*B_val*C_val
                Sp_site_ind_Sm_site_ind_p1(ind_i,ind_j) = result_ij
            
                ! For S-S+:
                if (b_i == 2 .AND. b_j == 1) then
                    B_val = 1.0d0
                else
                    B_val = 0.0d0
                end if
                result_ij = A_val*B_val*C_val
                Sm_site_ind_Sp_site_ind_p1(ind_i,ind_j) = result_ij
            
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
                result_ij = A_val*B_val*C_val
                Sz_site_ind_Sz_site_ind_p1(ind_i,ind_j) = result_ij
            
            end do
            end do
        
            H_full = H_full + 1.0d0/2.0d0 *(Sp_site_ind_Sm_site_ind_p1+Sm_site_ind_Sp_site_ind_p1) + Sz_site_ind_Sz_site_ind_p1
            
        end do
        
    
        write(*,*) 'H_full matrix diagonalization'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        !N_mat = N_spin_max ! The order of the matrix
        vl = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il=1  !the index of the smallest eigenvalue to be returned.
        iu=N_spin_max !the index of the largest eigenvalue to be returned.
        abstol=10**(-10.0d0)
        ldz=N_spin_max
        liwork= 10*N_spin_max
        lwork = 2*N_spin_max !check with work complex(1)
        lrwork= 24*N_spin_max

        allocate (eigen_vec(N_spin_max,N_spin_max), work_complex(lwork),iwork(liwork),w(N_spin_max),isuppz(2*N_spin_max), rwork(lrwork))
        write(*,*) 'just before allocation'
        call zheevr(jobz, range, uplo, N_spin_max, H_full, N_spin_max, vl, vu, il, iu, abstol, m_eig, w, eigen_vec, ldz, isuppz, work_complex, lwork, rwork, lrwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work_complex(1)
        write(*,*) "number of eigeval found:", m_eig
        do i= 1, m_eig
            write(*,*) i, 'th eigenval is', w(i)
        end do
    
        write(*,*) 'eigenvec norm test:'
        do i=1,m_eig
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1, N_spin_max
                norm = norm + eigen_vec(j,i)*conjg(eigen_vec(j,i))
            end do
            write(*,*) i, norm
        end do
        deallocate ( eigen_vec, work_complex, iwork, w, isuppz, rwork  )
    
        deallocate( H_full, Sp_site_ind_Sm_site_ind_p1, Sm_site_ind_Sp_site_ind_p1, Sz_site_ind_Sz_site_ind_p1 )
        
        close(10)


    end subroutine H_XXX_full

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
                           
        do site_ind = 1, (N_spin-1)
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

    
    subroutine H_XXX_diag(N_spin, J_spin)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none

        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer :: i, j, ind_i, ind_j, N_spin_max
        double precision, allocatable :: H_full(:,:)
        double precision :: norm, H_full_ij, J_spin_number
 
        CHARACTER*1 :: jobz, range, uplo
        integer :: il, iu, ldz, liwork, lwork, info, m_eig
        double precision :: vl, vu, abstol
        double precision, allocatable :: work(:), eigen_vec(:,:)
        integer, allocatable :: iwork(:), isuppz(:)
        double precision, allocatable :: w(:)
        
        CHARACTER(len=100) :: file_name1, file_name2, file_name3
        CHARACTER(len=10) :: N_spin_charachter

        ! Write the integer into a string:
        write(N_spin_charachter, '(i0)') N_spin

        file_name1 = 'H_eigenvals_full_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name2 = 'H_eigenvectors_full_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name3 = 'H_norm_test_ful_' // trim(adjustl(N_spin_charachter)) // '.dat'

        open (unit=10, file= trim(file_name1), recl=512)
        open (unit=11, file=trim(file_name2), recl=512)
        open (unit=12, file= trim(file_name3), recl=512)
        
        ! open (unit=11, file="H_full.dat", recl=512)
        ! lets rewrite our matlab code
        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
    
        call omp_mkl_small_test()
        
        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
        
        allocate( H_full(N_spin_max, N_spin_max))
        J_spin_number = J_spin
        H_full = 0.0d0

        !!! H_full filling
        write(*,*) 'H_full: '
          
        !$OMP PARALLEL DO
        do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max
                call H_XXX_filling(N_spin, J_spin_number, ind_i, ind_j, H_full_ij)
                H_full(ind_i,ind_j) = H_full_ij
                ! write(11,*) ind_i, ind_j, H_full(ind_i,ind_j)
            end do
        end do
        !$OMP END PARALLEL DO
        
    
        write(*,*) 'H_full matrix diagonalization'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        !N_mat = N_spin_max ! The order of the matrix
        vl     = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu     = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il     = 1  !the index of the smallest eigenvalue to be returned.
        iu     = N_spin_max !the index of the largest eigenvalue to be returned.
        abstol = 10**(-10.0d0)
        ldz    = N_spin_max
        lwork  = 26*N_spin_max
        liwork = 10*N_spin_max

        allocate ( w(N_spin_max), eigen_vec(N_spin_max,N_spin_max), isuppz(2*N_spin_max), work(lwork), iwork(liwork) )
        write(*,*) 'just before allocation'
       
        call dsyevr(jobz, range, uplo, N_spin_max, H_full, N_spin_max, vl, vu, il, iu, abstol, m_eig, w, eigen_vec, N_spin_max, isuppz, work, lwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work(1)
        write(*,*) "number of eigeval found:", m_eig

        write(10,*) 'i-th , eigenvalue'
        do i= 1, m_eig
            write(10,*) i, ',', w(i) 
        end do

        write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        do i=1, m_eig
            do j=1, N_spin_max
                write(11,*) eigen_vec(j,i)
            end do
            write(11,*) " "
        end do

        write(*,*) 'eigenvec norm test:'
        do i=1,m_eig
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1, N_spin_max
                norm = norm + eigen_vec(j,i)*(eigen_vec(j,i))
            end do
            write(12,*) i, norm
        end do
        
        deallocate (H_full,  w, eigen_vec, isuppz, work, iwork )
                            
        close(10)
        close(11)
        close(12)

    end subroutine H_XXX_diag

    subroutine H_XXX_diag_with_target_dense(N_spin, J_spin, hash)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none

        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer, allocatable, intent(in) :: hash(:)
        integer :: i, j, ind_i, ind_j, N_spin_max, max_val_hash
        double precision, allocatable :: H_full(:,:)
        double precision :: norm, H_full_ij, J_spin_number
 
        CHARACTER*1 :: jobz, range, uplo
        integer :: il, iu, ldz, liwork, lwork, info, m_eig
        double precision :: vl, vu, abstol
        double precision, allocatable :: work(:), eigen_vec(:,:)
        integer, allocatable :: iwork(:), isuppz(:)
        double precision, allocatable :: w(:)
        
        CHARACTER(len=100) :: file_name1, file_name2, file_name3
        CHARACTER(len=10) :: N_spin_charachter

        ! Write the integer into a string:
        write(N_spin_charachter, '(i0)') N_spin

        file_name1 = 'H_eigenvals_full_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name2 = 'H_eigenvectors_full_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name3 = 'H_norm_test_ful_' // trim(adjustl(N_spin_charachter)) // '.dat'

        open (unit=10, file= trim(file_name1), recl=512)
        open (unit=11, file=trim(file_name2), recl=512)
        open (unit=12, file= trim(file_name3), recl=512)
        
        ! open (unit=11, file="H_full.dat", recl=512)
        ! lets rewrite our matlab code
        ! N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
    
        !call omp_mkl_small_test()
        
        N_spin_max = 2**N_spin
        max_val_hash = maxval(hash)
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
        
        allocate( H_full(N_spin_max, N_spin_max))

        J_spin_number = J_spin
        H_full = 0.0d0

        !!! H_full filling
        write(*,*) 'H_full: '
          
        ! !$OMP PARALLEL DO
        do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max
                H_full_ij = 0.0d0
                
                if (hash(ind_i) > 0 .AND. hash(ind_j) > 0) then
                    call H_XXX_filling(N_spin, J_spin, hash(ind_i), hash(ind_j), H_full_ij)
                endif

                H_full(ind_i,ind_j) = H_full_ij
                ! write(11,*) ind_i, ind_j, H_full(ind_i,ind_j)
            end do
        end do
        ! !$OMP END PARALLEL DO

        write(*,*) 'H matrix for target 0 '
        do i = 1, N_spin_max
            do j = 1, N_spin_max

                write(*, '(I2,X)', advance='no'), H_full(i,j)

            enddo 
            write(*,*) ' ' ! this gives you the line break
        enddo
        
    
        write(*,*) 'H_full matrix diagonalization'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        !N_mat = N_spin_max ! The order of the matrix
        vl     = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu     = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il     = 1  !the index of the smallest eigenvalue to be returned.
        iu     = N_spin_max !the index of the largest eigenvalue to be returned.
        abstol = 10**(-10.0d0)
        ldz    = N_spin_max
        lwork  = 26*N_spin_max
        liwork = 10*N_spin_max

        allocate ( w(N_spin_max), eigen_vec(N_spin_max,N_spin_max), isuppz(2*N_spin_max), work(lwork), iwork(liwork) )
        write(*,*) 'just before allocation'
       
        call dsyevr(jobz, range, uplo, N_spin_max, H_full, N_spin_max, vl, vu, il, iu, abstol, m_eig, w, eigen_vec, N_spin_max, isuppz, work, lwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work(1)
        write(*,*) "number of eigeval found:", m_eig

        write(10,*) 'i-th , eigenvalue'
        do i= 1, m_eig
            write(10,*) i, ',', w(i) 
        end do

        write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        do i=1, m_eig
            do j=1, N_spin_max
                write(11,*) eigen_vec(j,i)
            end do
            write(11,*) " "
        end do

        write(*,*) 'eigenvec norm test:'
        do i=1,m_eig
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1, N_spin_max
                norm = norm + eigen_vec(j,i)*(eigen_vec(j,i))
            end do
            write(12,*) i, norm
        end do
        
        deallocate (H_full,  w, eigen_vec, isuppz, work, iwork )
                            
        close(10)
        close(11)
        close(12)

    end subroutine H_XXX_diag_with_target_dense

    subroutine H_XXX_feast_vec_fill(N_spin, J_spin, no_of_nonzero)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none
        
        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer, intent(out) :: no_of_nonzero
        integer :: i, j, ind_i, ind_j,  N_spin_max, counter, ind_ia            
        double precision :: norm, H_full_ij

        open (unit=20, file="ia_h.dat", recl=512)
        open (unit=21, file="ja_h.dat", recl=512)
        open (unit=22, file="values_array_h.dat", recl=512)
  
        ! N spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1

        !N_spin = 4
        !J_spin = 1
        N_spin_max = 2**N_spin
        write(*,*) 'be carefull about N_spin_max max size'
        write(*,*) 'test for N_spin_max', N_spin_max
        
        ! we start with saving ia, ja and val_arr into files (needs to be compared with ram only approach)  
        
        ! DODAC TU OMP! nie zrobisz przy zapisie do pliku
       
        !!! H_full filling
        write(*,*) 'expected counter: =', N_spin_max*(N_spin_max+1)/2
        counter = 0
        ind_ia = 0
        do ind_i = 1, N_spin_max ! do wymiaru podprzestrzeni
            ind_ia = ind_ia + 1
            write(20,*) ind_ia
            
            do ind_j = ind_i, N_spin_max ! do wymiaru podprzestrzeni
                call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)

                !if hash(ind_i) > 0 .AND. hash(ind_j) > 0 then
                    ! call H_XXX_filling(N_spin, J_spin, hash(ind_i), hash(ind_j), H_full_ij)
                ! end if

                if (ind_i  == ind_j) then
                    !write(20, *)
                    write(21, *) ind_j
                    write(22, *) H_full_ij
                    counter = counter + 1
                else if (H_full_ij .NE. 0.0d0) then
                    !write(20, *) 
                    ind_ia = ind_ia + 1
                    write(21, *) ind_j
                    write(22, *) H_full_ij
                    counter = counter + 1
                else
                    counter = counter + 1
                end if
                
            end do
            
        end do
        
        write(20,*) ind_ia+1
        
        write(*,*) 'counter test = ', counter 
        no_of_nonzero = ind_ia
        write(*,*) 'number of non-zero elements = ', no_of_nonzero
        
        close(20)
        close(21)
        close(22)

    end subroutine H_XXX_feast_vec_fill

    
    subroutine H_XXX_feast_vec_diag(N_spin, no_of_nonzero)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none
                    
        integer, intent(in) :: N_spin, no_of_nonzero
        integer :: i, j, N_spin_max, info, m_eig, n, loop, m0, fpm(128), ja_temp, ia_temp

        CHARACTER*1 :: uplo
        double precision, allocatable :: values_array(:), x(:,:), e(:), res(:)
        integer, allocatable :: ia(:), ja(:)
        double precision :: emin, emax, epsout, norm, val_arr_temp

        CHARACTER(len=100) :: file_name1, file_name2, file_name3
        CHARACTER(len=10) :: N_spin_charachter

        !Write the integer into a string:
        write(N_spin_charachter, '(i0)') N_spin

        file_name1 = 'H_eigenvals_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'
        !file_name2 = 'H_eigenvectors_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'
        file_name3 = 'H_norm_test_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'

        open (unit=30, file= trim(file_name1), recl=512)
        !open (unit=31, file=trim(file_name2), recl=512)
        open (unit=32, file= trim(file_name3) , recl=512)

        !FEAST diagonalization
        ! we need to store diagonal zeros as well !!!
        N_spin_max = 2**N_spin
        allocate( values_array(no_of_nonzero), ia(N_spin_max+1), ja(no_of_nonzero) )
        
        !reading of exsisting files
        open (unit=20, file="ia_h.dat", recl=512)
        open (unit=21, file="ja_h.dat", recl=512)
        open (unit=22, file="values_array_h.dat", recl=512)
        
        do i=1, N_spin_max+1           
            read(20, *) ia_temp
            ia(i) = ia_temp
        end do
            
        do i=1, no_of_nonzero
            read(22,*) val_arr_temp
            values_array(i) = val_arr_temp
            read(21, *) ja_temp
            ja(i) = ja_temp
        end do
        
        close(20)
        close(21)
        close(22)

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
        emin = -8.0d0 ! The lower ... &
        emax =  5.0d0  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = N_spin_max !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n. 
        !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]. 
        !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = N_spin_max
        allocate( x(n,m0), e(m0), res(m0) )
        write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
        write(*,*) 'Before dfeast_scsrev... '
        call dfeast_scsrev(uplo, n, values_array, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m_eig, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) ' dfeast_scsrev info=', info
        write(*,*) 'After  dfeast_scsrev... '
            
        if (info /= 0) then
        write(*,*) 'problem with  dfeast_scsrev, info=', info
        end if 
        
        write(*,*) ' dfeast_scsrev eigenvalues found= ', m_eig
        write(30,*) 'i-th , eigenvalue'
        do i = 1 , m_eig
            write(30,*) i, ',' , e(i)
        end do
        
        write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        do i=1, m_eig
            do j=1, n
                write(31,*) x(j,i)
            end do
            write(31,*) " "
        end do

        write(*,*) ' dfeast_scsrev eigenvec norm:'
        do i=1, m_eig
            norm = 0.0d0
            do j=1, n
                norm = norm + x(j,i)*(x(j,i))
            end do
            write(32,*) i, norm
        end do
                    
        deallocate(values_array, ia, ja, x, e, res)
        
        close(30)
        !close(31)
        close(32)

    end subroutine H_XXX_feast_vec_diag

    subroutine H_XXX_blocks_sz_feast()
        implicit none
        !basis_sz
        !HERE CREATE BLOCKS FOR H_XXX in that way to have them DIAGONALZIE VIA FEAST (in another subroutine)

    end subroutine H_XXX_blocks_sz_feast

    subroutine Rho_reduced_to_sz_basis
        implicit none

        !calculate density matrix rho from eigenvector due according to the truncated basis 

    end subroutine Rho_reduced_to_sz_basis

    subroutine Entropy_calculation()
        !size_reduced_basis, rho_reduced
        implicit none
        !calculate entropy according to reduced density matrix and normalazing it due to the truncated basis
        ! diagonalize reduced density matrix via LAPAC, check if its symmetric!

        !implicit none 
        !use omp_lib
        !use mkl_vsl
        !use omp_mkl_test

        !integer, intent(in) :: size_reduced_basis
        !double precision, allocatable, intent(in) :: rho_reduced

        !allocate(rho_reduced(size_reduced_basis))
 
        !def calculate_entropy(self,rho_reduced,n):
        
        !Here depending if s = 1/2 or s = 1 you need to change the base of log 
        !n - number of spins in the subsystem
        !eigen_rho, vectors = self.eig_diagonalize(rho_reduced) 
        
        !entropy = -sum(eigen_rho*np.log(eigen_rho, where=0<eigen_rho, out=0.0*eigen_rho))
        !eigen_rho_nonzero = eigen_rho[(eigen_rho > 10e-8) & (eigen_rho < 1.0)]
        !entropy = -np.sum(eigen_rho_nonzero * np.log2(eigen_rho_nonzero))
        
        !entropy = 0
        !for i in range(len(eigen_rho)):
            !print(eigen_rho[i])
        !    if eigen_rho[i] <= 10e-8:
        !        entropy += 0.0
                
         !   elif eigen_rho[i] == 1.0:
         !       entropy += 0.0
                
          !  else:
         !       entropy += -(eigen_rho[i]*np.log2(eigen_rho[i]))

        !deallocate(rho_reduced)
    
    end subroutine Entropy_calculation

end module spin_systems