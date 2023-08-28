include "mkl_vsl.f90" !necessary for intel mkl random numbers

module tests_module
    contains
    
    subroutine omp_mkl_small_test()
        use omp_lib
        use mkl_vsl
        implicit none

        integer :: i, j, omp_id, brng, rand_method, rand_seed, errcode, no_of_rand_num
        real(kind=8) rand_start, rand_end
        real(kind=8), allocatable :: rand_num(:)
        type (VSL_STREAM_STATE) :: stream(1)

        CHARACTER*1 :: jobz, range, uplo
        integer :: N_mat, il, iu, ldz, liwork, m, lwork, lrwork, info
        integer, allocatable :: iwork(:), isuppz(:)
        double precision :: vl, vu, abstol
        double complex, allocatable :: H_test(:,:), work_complex(:), eigen_vec(:,:)
        double precision, allocatable :: rwork(:), w(:)
        double complex :: norm

        double complex, allocatable  :: inv_H_test(:,:)
        integer, allocatable :: ipv(:)

        write(*,*) '----------- OMP MKL small matrix test:-----------'
        write(*,*) " "

        !   0 - part checking OpenMP + IntelMKL (oneAPI) random numbers/diag/inversion/
        write(*,*) 'OpenMP processors check'
        !$OMP PARALLEL PRIVATE(omp_id)
        omp_id = omp_get_thread_num()
        write(*,*) 'OpenMP thread check', omp_id
        !$OMP END PARALLEL

    

        write(*,*) 'Math Kernel Library random numbers check'
        ! *** random number generator parameters *** !
        brng   = VSL_BRNG_MT2203
        rand_method = VSL_RNG_METHOD_UNIFORM_STD
        rand_start = 0.0
        rand_end = 1.0
        rand_seed = 1
        no_of_rand_num = 10
        allocate(rand_num(no_of_rand_num))
        errcode = vslnewstream(stream(1), brng, rand_seed)
        write(*,*) 'random stream init code: ', errcode
        errcode=vdrnguniform( rand_method, stream(1), no_of_rand_num, rand_num, rand_start, rand_end )
        write(*,*) 'random numbers error code: ', errcode
        do i = 1, no_of_rand_num
            write(*,*) i, rand_num(i)
        end do
        deallocate(rand_num)
        errcode=vsldeletestream(stream(1))
        write(*,*) 'random stream kill code: ', errcode
        !end do
        write(*,*) 'Small Hermitean matrix diagonalization test'
        jobz  = 'V' !Compute eigenvalues and eigenvectors.
        range = 'I' !IL-th through IU-th eigenvalues will be found
        uplo  = 'U' !Upper triangle of A is stored;
        N_mat = 3 ! The order of the matrix
        vl = -30.0d0 ! lower bound of eigenvalues (for range='V')
        vu = 30.0d0  ! upper bound of eigenvalues (for range='V')
        il=1  !the index of the smallest eigenvalue to be returned.
        iu=N_mat !the index of the largest eigenvalue to be returned.
        abstol=10**(-10.0d0)
        ldz=N_mat
        liwork= 10*N_mat
        lwork=2*N_mat
        lrwork=24*N_mat
        allocate (H_test(N_mat,N_mat))
        H_test(1,1) = dcmplx(1.0d0, 0.0d0)
        H_test(1,2) = dcmplx(1.0d0, 1.0d0)
        H_test(1,3) = dcmplx(2.0d0, 2.0d0)
        H_test(2,1) = dcmplx(1.0d0,-1.0d0)
        H_test(2,2) = dcmplx(2.0d0, 0.0d0)
        H_test(2,3) = dcmplx(3.0d0, 3.0d0)
        H_test(3,1) = dcmplx(2.0d0,-2.0d0)
        H_test(3,2) = dcmplx(3.0d0,-3.0d0)
        H_test(3,3) = dcmplx(3.0d0, 0.0d0)
        allocate ( eigen_vec(N_mat,N_mat), work_complex(lwork),iwork(liwork),w(N_mat),isuppz(2*N_mat), rwork(lrwork)  )
        call zheevr(jobz, range, uplo, N_mat, H_test, N_mat, vl, vu, il, iu, abstol, m, w, eigen_vec, ldz, isuppz, work_complex, lwork, rwork, lrwork, iwork, liwork, info)
        write(*,*) "general zheevr info:", info
        write(*,*) lwork, " vs optimal lwork:", work_complex(1)
        write(*,*) "number of eigeval found:", m
        do i= 1, m
            write(*,*) i, 'th eigenval is', w(i)
        end do
        write(*,*) 'compare with Matlab: -2.4801476299237900'
        write(*,*) 'compare with Matlab:  0.5056116166572769'
        write(*,*) 'compare with Matlab:  7.9745360132665200'

    

        write(*,*) 'eigenvec norm test:'
        do i=1,m
            norm = dcmplx(0.0d0, 0.0d0)
            do j=1,N_mat
                norm = norm + eigen_vec(j,i)*conjg(eigen_vec(j,i))
            end do
            write(*,*) i, norm
        end do
        deallocate ( eigen_vec, work_complex, iwork, w, isuppz, rwork  )
        !
        write(*,*) 'Small Hermitean matrix inversion test'
        H_test(1,1) = dcmplx(1.0d0, 0.0d0)
        H_test(1,2) = dcmplx(1.0d0, 1.0d0)
        H_test(1,3) = dcmplx(2.0d0, 2.0d0)
        H_test(2,1) = dcmplx(1.0d0,-1.0d0)
        H_test(2,2) = dcmplx(2.0d0, 0.0d0)
        H_test(2,3) = dcmplx(3.0d0, 3.0d0)
        H_test(3,1) = dcmplx(2.0d0,-2.0d0)
        H_test(3,2) = dcmplx(3.0d0,-3.0d0)
        H_test(3,3) = dcmplx(3.0d0, 0.0d0)
        lwork = N_mat
        allocate(inv_H_test(N_mat,N_mat), ipv(N_mat), work_complex(lwork) )
        call zgetrf(N_mat, N_mat, H_test, N_mat, ipv, info ) !number of operations: (8/3)*N_mat^3
        write(*,*) 'blas LU info ', info
        call zgetri(N_mat, H_test, N_mat, ipv, work_complex, N_mat, info )
        !number of operations: (16/3)*M^3
        write(*,*) 'blas inversion info ', info
        write(*,*) lwork, 'inversion optimal lwork ', work_complex(1)
        do i = 1, N_mat
            do j = 1, N_mat
                write(*,*) i, j, H_test(i,j)
            end do
        end do
        write(*,*) 'inverted matrix should be:'
        write(*,*) ' 1.2+0.0i, -0.9+0.3i,  0.4-0.2i'
        write(*,*) '-0.9-0.3i,  0.5+0.0i, -0.1+0.3i'
        write(*,*) ' 0.4+0.2i, -0.1-0.3i,  0.0+0.0i'
        deallocate(H_test, ipv, work_complex)

        write(*,*) '----------- END OMP MKL test -----------'
        write(*,*) " "

    end subroutine omp_mkl_small_test

subroutine mmm_csr_test()

        use mkl_spblas
        use iso_C_binding
        implicit none 

        integer :: i, j, stat_permutation, stat_H_csr, stat_csr_mutliplication, info, stat_export, nrow, ncol, info_spmmd
        double precision, allocatable :: values_array_H(:), values_array_P(:)
        integer, allocatable :: ia_H_B(:), ia_H_E(:), ja_H(:), ia_P_B(:), ia_P_E(:), ja_P(:)

        !BLAS types
        type(sparse_matrix_t) :: p_matrix_csr
        type(sparse_matrix_t) :: H_matrix_csr
        type(sparse_matrix_t) :: H_matrix_permuted
        type(sparse_matrix_t) :: H_permuted
        type(matrix_descr) :: descr_H, descr_P
        !C_pointer types
        integer(C_INT) :: indexing
        type(C_PTR)   :: ia_start, ia_end, ja, values

        integer, pointer :: ia_export_B(:), ia_export_E(:), ja_export(:)
        double precision, pointer :: values_array_export(:)

        write(*,*) '----------- Sparse BLAS 4x4 matrix test:-----------'
        write(*,*) " "
        allocate(values_array_P(4), ia_P_B(4), ia_P_E(4), ja_P(4))

        ! P matrix
        ! 1 0 0 0
        ! 0 0 1 0
        ! 0 1 0 0
        ! 0 0 0 1 
        !values vector
        values_array_P(1) = 1.0d0
        values_array_P(2) = 1.0d0
        values_array_P(3) = 1.0d0
        values_array_P(4) = 1.0d0
        !columns vector
        ja_P(1) = 1
        ja_P(2) = 3
        ja_P(3) = 2
        ja_P(4) = 4
        ! csr pointer B
        ia_P_B(1) = 1;
        ia_P_B(2) = 2;
        ia_P_B(3) = 3;
        ia_P_B(4) = 4;
        ! csr pointer E
        ia_P_E(1) = 2;
        ia_P_E(2) = 3;
        ia_P_E(3) = 4;
        ia_P_E(4) = 5;


        ! H matrix symmetric upper triangler to csr3
        allocate(values_array_H(5), ia_H_B(4), ia_H_E(4), ja_H(5))
        ! H matrix
        ! 1 0 1 0
        ! 0 1 0 2
        ! 1 0 1 0
        ! 0 2 0 1
        !values vector
        values_array_H(1) = 1.0d0
        values_array_H(2) = 1.0d0
        values_array_H(3) = 1.0d0
        values_array_H(4) = 2.0d0
        values_array_H(5) = 1.0d0
        values_array_H(6) = 1.0d0
        !columns vector
        ja_H(1) = 1
        ja_H(2) = 3
        ja_H(3) = 2
        ja_H(4) = 4
        ja_H(5) = 3
        ja_H(6) = 4
        ! csr pointer B
        ia_H_B(1) = 1;
        ia_H_B(2) = 3;
        ia_H_B(3) = 5;
        ia_H_B(4) = 6;
        ! csr pointer E
        ia_H_E(1) = 3;
        ia_H_E(2) = 5;
        ia_H_E(3) = 6;
        ia_H_E(4) = 7;

        !rewriting permutation matrix to intel mkl internal csr format
        stat_permutation = mkl_sparse_d_create_csr(p_matrix_csr, SPARSE_INDEX_BASE_ONE, 4, 4 , ia_P_B, ia_P_E, ja_P, values_array_P)
        print *, "stat CSR Permutation matrix create = ", stat_permutation

        stat_export = mkl_sparse_d_export_csr(p_matrix_csr, indexing, nrow, ncol, ia_start, ia_end, ja, values)
        print *, "Export P = ", stat_export

        !   Converting C into Fortran pointers
        call C_F_POINTER(ia_start, ia_export_B, [nrow])
        call C_F_POINTER(ia_end  , ia_export_E  , [nrow]) 
        call C_F_POINTER(ja  , ja_export  , [ia_export_E(nrow)-indexing])
        call C_F_POINTER(values    , values_array_export    , [ia_export_E(nrow)-indexing])
        
        write(*,*)  ' '
        write(*,*) 'Exported CSR format of P: '
        print *, 'values array: ', values_array_export
        print *, 'ja_P: ', ja_export
        print *, 'ia_P; pointerB: ', ia_export_B
        print *, 'ia_P; pointerE: ', ia_export_E

        descr_P%type = SPARSE_MATRIX_TYPE_GENERAL
        info = MKL_SPARSE_OPTIMIZE(p_matrix_csr)
        write(*,*) 'info optimize P', info
    
        !rewriting H matrix to intel mkl internal csr format
        stat_H_csr = mkl_sparse_d_create_csr(H_matrix_csr, SPARSE_INDEX_BASE_ONE, 4, 4 , ia_H_B, ia_H_E, ja_H, values_array_H)
        print *, "stat CSR H matrix create = ", stat_H_csr

        descr_H % type = SPARSE_MATRIX_TYPE_SYMMETRIC ! it describes H matrix, P is general (full csr) by assumption of mkl_sparse_sypr
        descr_H % mode = SPARSE_FILL_MODE_UPPER
        descr_H % diag = SPARSE_DIAG_NON_UNIT

        info = MKL_SPARSE_OPTIMIZE(H_matrix_csr)
        write(*,*) 'info optimize H', info
        !info = mkl_sparse_order(H_matrix_csr)
        !write(*,*) 'info optimize P', info

        stat_export = mkl_sparse_d_export_csr(H_matrix_csr, indexing, nrow, ncol, ia_start, ia_end, ja, values)
        print *, "Export H_csr = ", stat_export

        !   Converting C into Fortran pointers
        call C_F_POINTER(ia_start, ia_export_B, [nrow])
        call C_F_POINTER(ia_end  , ia_export_E  , [nrow]) 
        call C_F_POINTER(ja  , ja_export  , [ia_export_E(nrow)-indexing])
        call C_F_POINTER(values    , values_array_export    , [ia_export_E(nrow)-indexing])
        
        write(*,*)  ' '
        write(*,*) 'Exported CSR format of H not permuted: '
        print *, 'values array: ', values_array_export
        print *, 'ja_H: ', ja_export
        print *, 'ia_H; pointerB: ', ia_export_B
        print *, 'ia_H; pointerE: ', ia_export_E


        !Doing a product of three csr matrices P @ H @ P.T <-SPARSE_OPERATION_NON_TRANSPOSE
        stat_csr_mutliplication = mkl_sparse_sypr(SPARSE_OPERATION_NON_TRANSPOSE, p_matrix_csr, H_matrix_csr, descr_H, H_matrix_permuted, SPARSE_STAGE_FULL_MULT)

        print *, "stat P @ H matrix @ P.T  = ", stat_csr_mutliplication

        info = mkl_sparse_order(H_matrix_permuted)
        stat_export = mkl_sparse_d_export_csr(H_matrix_permuted, indexing, nrow, ncol, ia_start, ia_end, ja, values)

        print *, "Export H permuted  = ", stat_export
        
        ! H permuted
        ! 1 1 0 0
        ! 1 1 0 0
        ! 0 0 1 2
        ! 0 0 2 1

        !   Converting C into Fortran pointers
        call C_F_POINTER(ia_start, ia_export_B, [nrow])
        call C_F_POINTER(ia_end  , ia_export_E  , [nrow]) 
        call C_F_POINTER(ja  , ja_export  , [ia_export_E(nrow)-indexing])
        call C_F_POINTER(values    , values_array_export    , [ia_export_E(nrow)-indexing])
        
        write(*,*)  ' '
        write(*,*) 'Exported CSR format of H permuted: '
        print *, 'values array: ', values_array_export
        print *, 'ja_H_permuted: ', ja_export
        print *, 'ia_H_permuted; pointerB: ', ia_export_B
        print *, 'ia_H_permuted; pointerE: ', ia_export_E
        
        !   Release internal representation of CSR matrix
        info = mkl_sparse_destroy(p_matrix_csr)
        info = mkl_sparse_destroy(H_matrix_csr)
        info = mkl_sparse_destroy(H_matrix_permuted)

        deallocate(values_array_P, ia_P_B, ia_P_E, ja_P)
        deallocate(values_array_H, ia_H_B, ia_H_E, ja_H)

        write(*,*) '----------- END Sparse BLAS test: -----------'
        write(*,*) " "

    end subroutine mmm_csr_test

    subroutine sparse_dfeast_test()
        
        implicit none
        double complex, allocatable :: test_matrix(:,:)
        CHARACTER*1 :: jobz, uplo            
        integer :: i, n, lda, lwork, info
        double precision, allocatable :: eigenvalues(:), rwork(:)
        double complex, allocatable ::  work(:,:)
        double precision, allocatable :: values_array(:), x(:,:)
        integer, allocatable :: ia(:), ja(:)
        integer :: fpm(128)
        double precision :: emin, emax
        integer :: m0
        double precision :: epsout
        integer :: loop
        double precision, allocatable :: e(:)
        integer :: m
        double precision, allocatable :: res(:)
        
        write(*,*) '----------- FEAST test START small matrix: -----------'
        write(*,*) " "

        n = 5
        allocate (test_matrix(n,n))
        test_matrix = dcmplx(0.0d0, 0.0d0)
        test_matrix(1,1) = dcmplx(1.0d0, 0.0d0)
        test_matrix(2,2) = dcmplx(0.0d0, 0.0d0)
        test_matrix(3,3) = dcmplx(2.0d0, 0.0d0)
        test_matrix(4,4) = dcmplx(0.0d0, 0.0d0)
        test_matrix(5,5) = dcmplx(3.0d0, 0.0d0)
        test_matrix(1,4) = dcmplx(5.0d0, 0.0d0)
        test_matrix(2,3) = dcmplx(6.0d0, 0.0d0)
        test_matrix(3,5) = dcmplx(7.0d0, 0.0d0)
        
        !standard diagonalziation
        jobz='V' !N - only eigenvalues, V - eigenvalues + eigenvectors
        uplo='U' !upper hermitean matrix part passed to routine
        n=5
        lda=n
        lwork=2*n-1  
        allocate(work(lwork,lwork),rwork(3*n-2),eigenvalues(n))
        
        call zheev(jobz, uplo, n, test_matrix, lda, eigenvalues, work, lwork, rwork, info)
        if (info /= 0) then
            write(*,*) 'problem with zheev, info=', info
        end if 
            
        do i=1,n
            write(*,*) i, eigenvalues(i) 
        enddo     
        deallocate(work,rwork,eigenvalues)

        write(*,*) 'Zheev diagonalization finished, FEAST start'

        !FEAST diagonalization
        ! we need to store diagonal zeros as well !!!
        allocate( values_array(8), ia(6), ja(8) )
         
        values_array(1) = dcmplx(1.0d0, 0.0d0)
        values_array(2) = dcmplx(5.0d0, 0.0d0)
        values_array(3) = dcmplx(0.0d0, 0.0d0)
        values_array(4) = dcmplx(6.0d0, 0.0d0)
        values_array(5) = dcmplx(2.0d0, 0.0d0)
        values_array(6) = dcmplx(7.0d0, 0.0d0)
        values_array(7) = dcmplx(0.0d0, 0.0d0)
        values_array(8) = dcmplx(3.0d0, 0.0d0)
        
        ja(1) = 1
        ja(2) = 4
        ja(3) = 2
        ja(4) = 3
        ja(5) = 3
        ja(6) = 5
        ja(7) = 4
        ja(8) = 5
        
        ia(1) = 1
        ia(2) = 3
        ia(3) = 5
        ia(4) = 7
        ia(5) = 8
        ia(6) = 9 !no of non-zero +1 !!!
           
        
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
        emin = -20.0d0 ! The lower ... &
        emax =  20.0d0  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = 5 !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n. 
           !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]. 
           !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = 5
        allocate( x(n,m0), e(m0), res(m0) )
        write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
        write(*,*) 'Before dfeast_hcsrev... '
        call dfeast_scsrev(uplo, n, values_array, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) 'dfeast_hcsrev info=', info
        write(*,*) 'After dfeast_hcsrev... '
            
        if (info /= 0) then
           write(*,*) 'problem with dfeast_hcsrev, info=', info
        end if 
        
        write(*,*) 'no of eigenvalues found= ', m
        do i = 1, m
            write(*,*) i, e(i)
        end do
                    
        deallocate(test_matrix, values_array, ia, ja, x, e, res)

        write(*,*) '----------- FEAST test END: -----------'
        write(*,*) " "

    end subroutine sparse_dfeast_test

    subroutine test_permutation_H_for_4_sites()
        use math_functions
        use spin_systems
        implicit none
        
        integer :: N_spin = 4
        double precision :: J_spin = 1.0
        integer :: i, j, N_spin_max
        integer, allocatable :: p_matrix(:,:), index_array_sz(:)
        double precision, allocatable :: H_full(:,:), H_permuted(:,:), Sp_site_ind_Sm_site_ind_p1(:,:), Sm_site_ind_Sp_site_ind_p1(:,:), Sz_site_ind_Sz_site_ind_p1(:,:)  
        integer :: site_ind, m, n, p, q, r, s
        integer :: ind_i, ind_j, a_i, a_j, b_i, b_j, c_i, c_j
        double precision :: H_full_ij

        write(*,*) '----------- TEST Hamiltonian permutation by matrix P defined by basis: -----------'
        write(*,*) "Test of permuting a H by P @ H @ P.T for N = 4 and J = 1.0"
        write(*,*) ' '

        N_spin_max = 2**N_spin

        allocate(index_array_sz(N_spin_max), H_full(N_spin_max, N_spin_max), H_permuted(N_spin_max, N_spin_max), Sp_site_ind_Sm_site_ind_p1(N_spin_max, N_spin_max), Sm_site_ind_Sp_site_ind_p1(N_spin_max, N_spin_max), Sz_site_ind_Sz_site_ind_p1(N_spin_max, N_spin_max), p_matrix(N_spin_max, N_spin_max))

        !$OMP PARALLEL DO
        do ind_i = 1, N_spin_max
            do ind_j = 1, N_spin_max
                call H_XXX_filling(N_spin, J_spin, ind_i, ind_j, H_full_ij)
                H_full(ind_i,ind_j) = H_full_ij
                ! write(11,*) ind_i, ind_j, H_full(ind_i,ind_j)
            end do
        end do
        !$OMP END PARALLEL DO

        !creating H for N = 4 
        call H_create_basis_sz(N_spin, index_array_sz)

        !creating P for N = 4 
        call Create_permutation_matrix(N_spin_max, index_array_sz, p_matrix)

        write(*,*) 'P matrix: '
        do i = 1, N_spin_max
            do j = 1, N_spin_max

                write(*, '(I2,X)', advance='no'), p_matrix(i,j)

            enddo 
            write(*,*) ' ' ! this gives you the line break
        enddo

        !permutating as P @ H @ P.T
        call mmm_mat_mul(H_permuted, dble(p_matrix), H_full, dble(p_matrix), 'n', 'n', 't')

        write(*,*) 'H matrix indicies for 4 sites: '
        do i = 1, N_spin_max
            do j = 1, N_spin_max

            print *,  i, j, H_permuted(i,j)

            enddo 
        enddo 

        write(*,*) 'H matrix for 4 sites: '
        do i = 1, N_spin_max
            do j = 1, N_spin_max

                write(*, '(I2,X)', advance='no'), H_permuted(i,j)

            enddo 
            write(*,*) ' ' ! this gives you the line break
        enddo

        write(*,*) '----------- END test Hamiltonian permutation by matrix P defined by basis: -----------'
        write(*,*)  " "
    
        deallocate(H_full, Sp_site_ind_Sm_site_ind_p1, Sm_site_ind_Sp_site_ind_p1, Sz_site_ind_Sz_site_ind_p1, p_matrix, H_permuted,index_array_sz)


    end subroutine test_permutation_H_for_4_sites

end module tests_module