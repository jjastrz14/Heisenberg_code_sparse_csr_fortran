module spin_systems
    contains

    subroutine H_create_basis_sz_sorted(N_spin, indices_Sz_basis_sorted)
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
        !commenting these lines
        !go to 90
        ! CHUNK OF CODE YOU WANT TO COMMENT OUT
        
        !sorting
        id = 'D'
        !dlasrt
        call dlapst(id, N_spin_max, Sz_basis, indices_Sz_basis_sorted, info ) !scalapack quicksort

        !subroutine sorts basis_sz array and creates array for permutation matrix
        ! call indexArrayReal(N_spin_max, basis_sz, index_array)
        !90 continue

        write(*,*) 'Indicies of Sorted Sz basis: '
        write(*,*) indices_Sz_basis_sorted

        write(*,*) 'Sorted Sz basis: '
        write(*,*) Sz_basis(indices_Sz_basis_sorted)

        deallocate(Sz_basis)

    end subroutine H_create_basis_sz_sorted


    subroutine H_create_basis_sz(N_spin, Sz_basis, basis_vector)
        !> Returns a binary representation of basis(i) as a string with at least N_spin bits. 
        !> abs(d) < 2^52
        ! dec2bin(d, n)
        use math_functions
        implicit none

        integer, intent(in) :: N_spin
        double precision :: Sz
        character(len=53), allocatable :: dec2bin(:)
        double precision, allocatable   :: basis_sz(:)
        double precision, allocatable, intent(out) :: Sz_basis(:)
        integer :: N_spin_max, i, j, info, N, interator
        integer, allocatable, intent(out) :: basis_vector(:,:)
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
        allocate (basis_vector(N_spin_max, N_spin))

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
        ! 0 == false => Sz=Sz+1/2, 1= true> = Sz=Sz-1/2 : spin_up -> 0, spin_down -> 1
        do i = 0 , N_spin_max-1
            basis(i) = i !here we want Sz from this i
            ! write(*,*) 'basis vector i =', i+1
            Sz = 0.0d0
            do j=0, N_spin-1
                bool = btest(i, j) 
                !write(*,*) j, bool

                if (bool==.FALSE.) then
                    Sz=Sz+1.0d0/2.0d0
                    ! this basis_vector is used for truncating basis for reduced density matrices 
                    basis_vector(i + 1, N_spin - j) = 0
                elseif (bool==.TRUE.) then
                    Sz=Sz-1.0d0/2.0d0
                    ! this basis_vector is used for truncating basis for reduced density matrices 
                    basis_vector(i + 1, N_spin - j) = 1
                else
                    write(*,*) 'something during basis creation wrong'
                    error stop
                end if  
            end do    
            Sz_basis(i+1) = Sz
            ! write(*,*) 'for basis vector i =', i, 'Sz = ', Sz_basis(i+1)
        end do

        ! this basis_vector is used for truncating basis for reduced density matrices 
        !write(*,*) 'Basis vector 0 and 1: ' 
        !do i = 1, N_spin_max 
          !  do j = 1, N_spin 
          !      write(*,*), i, j, basis_vector(i,j)
          !  end do 
        !end do 

        write(*,*) 'Summary of Sz basis: '
        do i=1,N_spin_max
            write(*,*) Sz_basis(i)
        end do

        deallocate(basis)

    end subroutine H_create_basis_sz

    subroutine Hash_basis_with_target(N_spin, target_sz, Sz_basis, basis_vector, hash, basis_rho_target)
        implicit none 

        integer, intent(in) :: target_sz, N_spin, basis_vector(:,:)
        double precision, intent(in) :: Sz_basis(:)
        integer, allocatable, intent(out) :: hash(:)
        integer, allocatable, intent(out) :: basis_rho_target(:,:)
        integer :: N, i,j, N_spin_max, max_val_hash_loc, min_val_hash_loc, max_val_hash

        N_spin_max = 2**N_spin
        allocate(hash(N_spin_max))
        
        N = 1
        do i = 1, N_spin_max
            if (Sz_basis(i) == target_sz) then
                !Sz_basis_target(i) = i
                !Sz_basis_target(N) = i
                hash(i) = N  
                N = N + 1
            else 
                hash(i) = -1 
            endif 
        end do 

        !write(*,*) 'Hash:'
        !write(*,*) hash

        max_val_hash_loc = MAXLOC(hash, dim = 1)
        min_val_hash_loc = FINDLOC(hash, 1,  dim = 1)
        max_val_hash = MAXVAl(hash)

        allocate(basis_rho_target(max_val_hash, N_spin))

        do i = 1, N_spin_max 
                if (hash(i) > 0) then
                    basis_rho_target(hash(i),:) = basis_vector(i,:)
                else if (hash(i) <= 0) then 
                    continue 
                else 
                    write(*,*) 'something wrong during hashing reduced basis procedure'
                    error stop
                end if 
        end do 

        !write(*,*) 'Hashed Basis vector 0 and 1: ' 
        !do i = 1, max_val_hash 
         !  do j = 1, N_spin 
         !     write(*,*), i, j, basis_rho_target(i,j)
         !  end do 
        !end do
        
        write(*,*) "Size of basis rho target", size(basis_rho_target, 1)

    end subroutine Hash_basis_with_target

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

        write(*,*) ' dsyevr eigenvectors to file :'
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

    subroutine H_XXX_block_diag_with_target_dense(N_spin, J_spin, hash)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none

        integer, intent(in) :: N_spin
        double precision, intent(in) :: J_spin
        integer, allocatable, intent(in) :: hash(:)
        integer :: i, j, ind_i, ind_j, N_spin_max, max_val_hash_loc, min_val_hash_loc, max_val_hash, dim
        double precision, allocatable :: H_full_block(:,:)
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
        ! Example for N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        
        N_spin_max = 2**N_spin
        max_val_hash_loc = MAXLOC(hash, dim = 1)
        min_val_hash_loc = FINDLOC(hash, 1,  dim = 1)
        max_val_hash = MAXVAl(hash)

        write(*,*) 'max val hash', max_val_hash
        write(*,*) 'max val hash location', max_val_hash_loc
        write(*,*) 'min val hash location', min_val_hash_loc
        
        allocate( H_full_block(max_val_hash, max_val_hash))

        J_spin_number = J_spin
        H_full_block = 0.0d0

        !!! H_full filling
        write(*,*) 'H block: '
          
        !$OMP PARALLEL DO
        do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max

                if (hash(ind_i) > 0 .AND. hash(ind_j) > 0) then
                    call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
                    H_full_block(hash(ind_i),hash(ind_j))= H_full_ij
                endif

            end do
        end do
        !$OMP END PARALLEL DO
        
        !printing above result

        do ind_i = 1, max_val_hash
            do ind_j = 1, max_val_hash

                print *,  ind_i, ind_j, H_full_block(ind_i,ind_j)

            end do
        end do
        

        write(*,*) 'H matrix for target 0 '
        do i = 1, max_val_hash
            do j = 1, max_val_hash

                write(*, '(I2,X)', advance='no'), H_full_block(i,j)

            enddo 
            write(*,*) ' ' ! this gives you the line break
        enddo
        
        
    
        write(*,*) 'H block matrix diagonalization'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        !N_mat = N_spin_max ! The order of the matrix
        vl     = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu     = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il     = 1  !the index of the smallest eigenvalue to be returned.
        iu     = max_val_hash !the index of the largest eigenvalue to be returned.
        abstol = 10**(-10.0d0)
        ldz    = N_spin_max
        lwork  = 26*N_spin_max
        liwork = 10*N_spin_max

        allocate ( w(max_val_hash), eigen_vec(max_val_hash,max_val_hash), isuppz(2*max_val_hash), work(lwork), iwork(liwork) )
        write(*,*) 'just before allocation'
       
        call dsyevr(jobz, range, uplo, max_val_hash, H_full_block, max_val_hash, vl, vu, il, iu, abstol, m_eig, w, eigen_vec, max_val_hash, isuppz, work, lwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work(1)
        write(*,*) "number of eigeval found:", m_eig

        write(10,*) 'i-th , eigenvalue'
        do i= 1, m_eig
            write(10,*) i, ',', w(i) 
            write(*,*) i, ',', w(i) 
        end do

        write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        do i=1, m_eig
            do j=1, max_val_hash
                write(11,*) eigen_vec(j,i)
            end do
            write(11,*) " "
        end do

        write(*,*) 'eigenvec norm test:'
        do i=1,m_eig
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1, max_val_hash
                norm = norm + eigen_vec(j,i)*(eigen_vec(j,i))
            end do
            write(12,*) i, norm
        end do
        
        deallocate (H_full_block, w, eigen_vec, isuppz, work, iwork )
                            
        close(10)
        close(11)
        close(12)

    end subroutine H_XXX_block_diag_with_target_dense

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

        !!! H_full filling
        write(*,*) 'expected counter: =', N_spin_max*(N_spin_max+1)/2
        counter = 0
        ind_ia = 0
        do ind_i = 1, N_spin_max ! do wymiaru podprzestrzeni
            ind_ia = ind_ia + 1
            write(20,*) ind_ia
            
            do ind_j = ind_i, N_spin_max ! do wymiaru podprzestrzeni
                call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)

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


    subroutine H_XXX_block_feast_vec_fill(N_spin, J_spin, hash, no_of_nonzero)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none
        
        integer, intent(in) :: N_spin
        integer, allocatable, intent(in) :: hash(:)
        double precision, intent(in) :: J_spin
        integer, intent(out) :: no_of_nonzero
        integer :: i, j, ind_i, ind_j,  N_spin_max, counter, ind_ia, max_val_hash_loc , min_val_hash_loc, max_val_hash    
        double precision :: norm, H_block_ij
        CHARACTER(len=10) :: N_spin_charachter

        write(N_spin_charachter, '(i0)') N_spin

        open (unit=20, file='ia_h_block' // trim(adjustl(N_spin_charachter)) // '.dat', recl=512)
        open (unit=21, file='ja_h_block' // trim(adjustl(N_spin_charachter)) // '.dat', recl=512)
        open (unit=22, file='values_array_h_block' // trim(adjustl(N_spin_charachter)) // '.dat', recl=512)
  
        ! N spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1

         ! Example for N=4 spin-1/2 Heisenberg XXX (Jx=Jy=Jz=J) model, J=1
        
        N_spin_max = 2**N_spin
        max_val_hash_loc = MAXLOC(hash, dim = 1)
        min_val_hash_loc = FINDLOC(hash, 1,  dim = 1)
        max_val_hash = MAXVAl(hash)

        write(*,*) 'max val hash', max_val_hash
        !write(*,*) 'max val hash location', max_val_hash_loc
        !write(*,*) 'min val hash location', min_val_hash_loc

        !!! H_full filling
        !write(*,*) 'H block feast vector fillings for target :'
          
        ! !$OMP PARALLEL DO
        !do ind_i = 1, N_spin_max
         !   do ind_j = 1, N_spin_max

          !!      if (hash(ind_i) > 0 .AND. hash(ind_j) > 0) then
          !          call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_block_ij)
          !          H_full_block(hash(ind_i),hash(ind_j))= H_block_ij
           !     endif
    !
           ! end do
        !end do
       ! !$OMP END PARALLEL DO

        write(*,*) 'expected counter: =', max_val_hash*(max_val_hash+1)/2
        counter = 0
        ind_ia = 0

        ! !$OMP PARALLEL DO

        do ind_i = 1, N_spin_max ! do wymiaru podprzestrzeni

            if (hash(ind_i) > 0) then
                ind_ia = ind_ia + 1
                write(20,*) ind_ia
            end if 

            do ind_j = ind_i, N_spin_max ! do wymiaru podprzestrzeni

                if (hash(ind_i) > 0 .AND. hash(ind_j) > 0) then
                    call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_block_ij)

                    if (ind_i  == ind_j) then
                        !write(20, *)
                        write(21, *) hash(ind_j)
                        write(22, *) H_block_ij
                        counter = counter + 1
                    else if (H_block_ij .NE. 0.0d0) then
                        !write(20, *) 
                        ind_ia = ind_ia + 1
                        write(21, *) hash(ind_j)
                        write(22, *) H_block_ij
                        counter = counter + 1
                    else 
                        counter = counter + 1
                    end if

                endif
            end do
        end do

        ! !$OMP END PARALLEL DO
        
        write(20,*) ind_ia+1
        
        write(*,*) 'counter test = ', counter 
        no_of_nonzero = ind_ia
        write(*,*) 'number of non-zero elements = ', no_of_nonzero
        
        close(20)
        close(21)
        close(22)

    end subroutine H_XXX_block_feast_vec_fill
 
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

    subroutine H_XXX_block_feast_vec_diag(N_spin, no_of_nonzero, hash, e, x)
        use omp_lib
        use mkl_vsl
        use tests_module
        implicit none
                    
        integer, intent(in) :: N_spin, no_of_nonzero
        integer, allocatable, intent(in) :: hash(:)
        integer :: i, j, N_spin_max, info, m_eig, n, loop, m0, fpm(128), ja_temp, ia_temp, max_val_hash
        CHARACTER*1 :: uplo
        double precision, allocatable :: values_array(:), res(:)
        double precision, allocatable, intent(out) :: x(:,:), e(:)
        integer, allocatable :: ia(:), ja(:)
        double precision :: emin, emax, epsout, norm, val_arr_temp

        CHARACTER(len=100) :: file_name1, file_name2, file_name3
        CHARACTER(len=10) :: N_spin_charachter

        !Write the integer into a string:
        write(N_spin_charachter, '(i0)') N_spin

       ! file_name1 = 'H_eigenvals_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'
        !file_name2 = 'H_eigenvectors_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'
        !file_name3 = 'H_norm_test_feast_' // trim(adjustl(N_spin_charachter)) // '.dat'

        !open (unit=30, file= trim(file_name1), recl=512)
        !open (unit=31, file=trim(file_name2), recl=512)
        !open (unit=32, file= trim(file_name3) , recl=512)

        !FEAST diagonalization
        ! we need to store diagonal zeros as well !!!
        N_spin_max = 2**N_spin
        max_val_hash = MAXVAl(hash)
        allocate( values_array(no_of_nonzero), ia(max_val_hash+1), ja(no_of_nonzero) )
        
        !reading of exsisting files
        open (unit=20, file='ia_h_block' // trim(adjustl(N_spin_charachter)) // '.dat', recl=512)
        open (unit=21, file='ja_h_block' // trim(adjustl(N_spin_charachter)) // '.dat', recl=512)
        open (unit=22, file='values_array_h_block' // trim(adjustl(N_spin_charachter)) // '.dat', recl=512)
  
        
        do i=1, max_val_hash+1           
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
        emin = -15d0 ! The lower ... &
        emax =  6d0  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = max_val_hash !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n. 
        !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]. 
        !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = max_val_hash
        allocate( x(n,m0), e(m0), res(m0) )
        !write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
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
        !write(30,*) 'i-th , eigenvalue'
        !do i = 1 , m_eig
           ! write(30,*) i, ',' , e(i)
        !end do
        
        !write(*,*) ' dfeast_scsrev eigenvectors to file :'
        ! saving of eigenvectors to file
        !do i=1, m_eig
           ! do j=1, n
            !    write(31,*) x(j,i)
            !end do
            !write(31,*) " "
        !end do

        !write(*,*) ' dfeast_scsrev eigenvec norm:'
        !do i=1, m_eig
            !norm = 0.0d0
            !do j=1, n
             !   norm = norm + x(j,i)*(x(j,i))
            !end do
            !write(32,*) i, norm
        !end do
                    
        deallocate(values_array, ia, ja, res)
        
        !close(30)
        !close(31)
        !close(32)

    end subroutine H_XXX_block_feast_vec_diag

    subroutine Rho_reduced_calculation(N_spin, basis_rho_target, size_of_sub_A, size_of_sub_B, eigen_vectors, index_energy, rho_reduced)
        use math_functions
        implicit none
        !calculate density matrix rho from eigenvector due according to the truncated basis 
          
        integer, intent(in) :: size_of_sub_A, size_of_sub_B, N_spin
        integer, allocatable, intent(in) :: basis_rho_target(:,:)
        double precision, dimension(:,:), allocatable, intent(in):: eigen_vectors
        integer, intent(in) :: index_energy
        double precision, dimension(:,:), allocatable, intent(out) :: rho_reduced

        double precision, dimension(:,:), allocatable :: psi
        integer, allocatable :: subsystem_A(:,:), subsystem_B(:,:), subsystem_A_set(:,:), subsystem_B_set(:,:)
        integer, dimension(2) :: new_entry
        integer, dimension(:,:), allocatable :: new_basis_rho_reduced
        integer, dimension(:,:), allocatable :: k_A, k_B
        integer :: i, j, k, l, index_new_basis_i, index_new_basis_j
        double precision :: v, trace
        logical :: found = .false.
        logical :: is_present
        
        ! here need to introduce this when you want divide subsystem not equally
        ! basis_combination_A = 2**size_of_sub_A
        ! basis_combination_B = 2**size_of_sub_B

        allocate(subsystem_A(size(basis_rho_target, 1) , size_of_sub_A), subsystem_B(size(basis_rho_target, 1) , size_of_sub_B))
        allocate(new_basis_rho_reduced(size(basis_rho_target, 1), 2)) !only i and j indicies

        ! Calculate the bases of subsystems A and B
        !print *, "Size of basis rho target", size(basis_rho_target, 1)

        do i = 1, size(basis_rho_target, 1)
            subsystem_A(i , :) = basis_rho_target(i , 1:size_of_sub_A)
            subsystem_B(i , :) = basis_rho_target(i , size_of_sub_B+1:N_spin)
        end do

        call RemoveDuplicates(subsystem_A, subsystem_A_set)
        call RemoveDuplicates(subsystem_B, subsystem_B_set)
    
        !print*, size(subsystem_A_set, 1)
        !print*, size(subsystem_B_set, 1)

        !subsystem_A_set(1,1) = 0
        !subsystem_A_set(1,2) = 1

        !subsystem_A_set(2,1) = 0
        !subsystem_A_set(2,2) = 0

        !subsystem_A_set(3,1) = 1
        !subsystem_A_set(3,2) = 0

        !subsystem_B_set(1,1) = 0
        !subsystem_B_set(1,2) = 0

        !subsystem_B_set(2,1) = 0
        !subsystem_B_set(2,2) = 1

        !subsystem_B_set(3,1) = 1
        !subsystem_B_set(3,2) = 0
        

        do k = 1, size(basis_rho_target, 1)
            !print*, subsystem_A(k , :)
            !print*, subsystem_B(k , :)

            call FindRowIndex(subsystem_A(k , :), subsystem_A_set, index_new_basis_i) !working
            call FindRowIndex(subsystem_B(k , :), subsystem_B_set, index_new_basis_j)

            !print *, " i ",  index_new_basis_i
            !print *, " j ", index_new_basis_j

           ! Calculate indicies for the new basis of rho reduced
            new_basis_rho_reduced(k,1) = index_new_basis_i
            new_basis_rho_reduced(k,2) = index_new_basis_j

        end do 

        ! print bases for subsystems 
        !print *, "Basis for whole rho:"
        !do i = 1, size(basis_rho_target, 1) 
         !   do j = 1, N_spin
         !       write(*,*), i, j, basis_rho_target(i,j)
         !   end do 
         !   print *, " "
        !end do 

        !print *, "Basis for subsystem A:"
        !do i = 1, size(basis_rho_target, 1)
        !    do j = 1, size_of_sub_A
        !       write(*,*), i, j, subsystem_A(i,j)
         !   end do 
        !    print *, " "
        !end do 

       ! print *, "SET subsystem A :"
        !do i = 1, size(subsystem_A_set, 1)
        !    do j = 1, size_of_sub_A
         !      write(*,*), i, j, subsystem_A_set(i,j)
        !    end do 
        !    print *, " "
        !end do 
        
        !print *, "Basis for subsystem B:"
        !do i = 1, size(basis_rho_target, 1)
         !   do j = 1, size_of_sub_B
         !       write(*,*), i, j, subsystem_B(i,j)
         !   end do 
        !   print *, " "
        !end do
        
       ! print *, "SET subsystem B :"
        !do i = 1, size(subsystem_B_set, 1)
        !    do j = 1, size_of_sub_B
        !       write(*,*), i, j, subsystem_B_set(i,j)
        !    end do 
        !   print *, " "
        !end do 

        !print *, "New basis indices : "
        !do i = 1, size(new_basis_rho_reduced, 1)
        !    do j = 1, 2
         !       write(*,*), i, j, new_basis_rho_reduced(i,j)
         !   end do 
         !   print *, " "
        !end do 

        deallocate(subsystem_A, subsystem_B)

        !create psi 
        allocate(psi(size(subsystem_A_set,1), size(subsystem_B_set,1)))
        psi = 0.0d0

        !print*, eigen_vectors(:,index_energy)

        do k = 1, size(eigen_vectors(:,index_energy))
        

            v = eigen_vectors(k,index_energy)
            !print *, "This is value of psi ", v
            psi(new_basis_rho_reduced(k,1) , new_basis_rho_reduced(k, 2)) = v

        end do 

        ! print*, "this is psi for index_energy", index_energy
        ! do i = 1, size(psi, 1)
        !   do j = 1, size(psi, 2)
        !       write(*,*), i, j, psi(i,j)
        !   end do 
        !end do 

        allocate(rho_reduced(size(psi,1),size(psi,1)))
        rho_reduced = 0.0d0
        
        ! for complex psi it should be also psi.conj 
        ! numpy dot product of 2 matrices from python is a matmul 
        !print *, "Rho reduced calculation: "
        ! this is |psi> <psi| = rho
        rho_reduced = MATMUL(psi, transpose(psi))

        !trace of rho_reduced check
        trace = 0.0d0
        do i = 1, size(rho_reduced, 1)
            trace = trace + rho_reduced(i,i)
        end do

        !do i = 1, size(rho_reduced, 1)
         !   do j = 1, size(rho_reduced, 2)
         !       write(*,*), i, j, rho_reduced(i,j)
         !   end do 
      ! end do 

        if (.not. (.999999999d0 <= trace .and. trace <= 1.000000001d0)) then
            print*, "Trace of the reduced density matrix is not equal to 1 ! : ", trace 
        end if 

    end subroutine Rho_reduced_calculation

    subroutine Entropy_calculation(systemA_size, size_reduced_basis, rho_reduced, entropy_value, eigen_value_check)
        !size_reduced_basis, rho_reduced
        implicit none
        !calculate entropy according to reduced density matrix and normalazing it due to the truncated basis
        ! diagonalize reduced density matrix using LAPACK, check if its symmetric!

        integer, intent(in) :: size_reduced_basis, systemA_size
        double precision, allocatable, intent(in) :: rho_reduced(:,:)
        double precision, intent(out):: entropy_value, eigen_value_check
        integer :: ind_i, i

        !diagonalization of non-symmetric reduced rho (reduced density matrix) using dgees from intem MKL
        CHARACTER*1 :: jobvs, sort
        integer :: il, iu, lda,ldvs, liwork, lwork, info, m_eig, n
        double precision :: vl, vu, abstol
        double precision, allocatable :: work(:), eigen_vec(:,:), vs(:)
        integer, allocatable :: iwork(:), isuppz(:)
        double precision, allocatable :: eigen_values_r(:), eigen_values_c(:)
        logical :: bwork, select

        !write(*,*) 'Rho_reduced matrix diagonalization'
        !write(*,*) 'check for bigger N if lower and upper bound for eigenvalues is correctly set'
        jobvs = 'N' !Compute eigenvalues and eigenvectors, NOT compute Shur vectors
        sort = 'N' !Do not sort eigenvalues
        n = size_reduced_basis ! The order of the matrix
        !vl     = -30.0d0 ! lower bound of eigenvalues (for range='V')
        !vu     = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il     = 1  !the index of the smallest eigenvalue to be returned.
        iu     = size_reduced_basis !the index of the largest eigenvalue to be returned.
        abstol = 10**(-10.0d0)
        lda    = size_reduced_basis
        ldvs   = size_reduced_basis
        lwork  = 10*size_reduced_basis
        liwork = 10*size_reduced_basis

        allocate ( eigen_values_r(size_reduced_basis), eigen_values_c(size_reduced_basis), eigen_vec(size_reduced_basis,size_reduced_basis), work(lwork), iwork(liwork) )
        !write(*,*) 'just before allocation'
        ! NOTICE WRONG PROCEDURE: RHEDUCED RHO ARE NON SYMMETRIC SO USE whole array diagonalization method
        ! now dgees is good
        call dgees(jobvs, sort, select, n, rho_reduced, lda, m_eig, eigen_values_r, eigen_values_c, vs, ldvs, work, lwork, bwork, info)

        !call gees(rho_reduced, eigen_values_r, eigen_values_c, [,vs] [,select] [,sdim] [,info])

        !dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
        !call dsyevr(jobz, range, uplo, size_reduced_basis, rho_reduced, size_reduced_basis, vl, vu, il, iu, abstol, m_eig, eigen_values, eigen_vec, size_reduced_basis, isuppz, work, lwork, iwork, liwork, info)
        !write(*,*) "general zheevr info:", info
        !write(*,*) lwork, " vs optimal lwork:", work(1)
        !write(*,*) "number of eigeval found in rho_reduced diag:", size(eigen_values_r, 1)

        !write(*,*) "eigenvalues of rho", eigen_values_r
        !write(*,*) " This is real eigenvalues part ", eigen_values_r

        do i = 1, m_eig
            if (.not. (.00000000d0 <= eigen_values_c(i) .and. eigen_values_c(i) <= 0.000000000001d0)) then
                print*, "Non zero complex eigenvalue in rho_reduced ", eigen_values_c(i)
            end if 
        end do 

        ! Entropy calculation
        entropy_value = 0.0d0
        eigen_value_check = 0.0d0

        do ind_i = 1, size(eigen_values_r, 1)
            if (eigen_values_r(ind_i) <= 10E-8) then 
                entropy_value = entropy_value + 0.0d0
            else if (eigen_values_r(ind_i) == 1.0d0 ) then
                entropy_value = entropy_value + 0.0d0
            else 
                entropy_value = entropy_value -(eigen_values_r(ind_i) * (log(eigen_values_r(ind_i))))
                eigen_value_check =  eigen_value_check + eigen_values_r(ind_i)
            end if 
        end do 

        !entropy_value = entropy_value
        entropy_value = entropy_value/(dble(systemA_size)*log(dble(systemA_size)))

        !write(*,*), "This is lambdas value check: ", eigen_value_check
        !write(*,*) "This is entropy for this energy", entropy_value
        !write(*,*) "This is entropy for this energy normalized", entropy_value/(dble(systemA_size)*log(dble(systemA_size))/log(2.0))
    
    end subroutine Entropy_calculation

end module spin_systems