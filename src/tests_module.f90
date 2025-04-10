include "mkl_vsl.f90" !necessary for intel mkl random numbers

module tests_module
    implicit none
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


    subroutine sparse_zfeast_test()
    
        implicit none
        double complex, allocatable :: test_matrix(:,:)
        CHARACTER*1 :: jobz, uplo            
        integer :: i, n, lda, lwork, info
        double complex, allocatable :: rwork(:)
        double complex, allocatable ::  work(:,:)
        double complex, allocatable :: values_array(:)
        integer, allocatable :: ia(:), ja(:)
        integer, dimension(64) :: fpm
        double precision :: emin, emax
        integer :: m0
        double precision :: epsout
        integer :: loop
        double precision, allocatable :: e(:), eigenvalues(:), x(:,:)
        integer :: m
        double precision, allocatable :: res(:)
        
        write(*,*) '----------- Complex FEAST test START small matrix: -----------'
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
        call feastinit(fpm)  ! function specifying default parameters fpm of FEAST algorithm
        fpm(1) = 1 !if 1 print runtime comments
        !fpm(2) = 12 !can be more, can be less
        !fpm(3) = 8 !eps 10^-fpm(3)
        !fpm(4) = 20 ! max number of feast loops
        !fpm(5) = 0 !initial subspace
        !fpm(6) = 0! stopping criterion
        !fpm(7) = 5 !Error trace sigle prec stop crit
        !fpm(14) = 0! standard use of feast
        !fpm(27) = 1 !check input matrices
        !fpm(28) = 1 !check if B is positive definite?    
        !fpm(13) = 0     
        uplo='U' ! If uplo = 'U', a stores the upper triangular parts of A.
        emin = -20.0d0 ! The lower ... &
        emax =  20.0d0  !  and upper bounds of the interval to be searched for eigenvalues
        m0 = 5 !On entry, specifies the initial guess for subspace dimension to be used, 0 < m0≤n. 
        !Set m0 ≥ m where m is the total number of eigenvalues located in the interval [emin, emax]. 
        !If the initial guess is wrong, Extended Eigensolver routines return info=3.
        n = 5
        allocate( x(n,m0), e(m0), res(m0) )
        !write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'

        !zfeast_hcsrev(UPLO,N,A,IA,JA,fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info)

        write(*,*) 'Before zfeast_hcsrev... '
        call zfeast_hcsrev(uplo, n, values_array, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) 'zfeast_hcsrev info=', info
        write(*,*) 'After zfeast_hcsrev... '
            
        if (info /= 0) then
        write(*,*) 'problem with zfeast_hcsrev, info=', info
        end if 
        
        write(*,*) 'no of eigenvalues found= ', m
        do i = 1, m
            write(*,*) i, e(i)
        end do
                    
        deallocate(test_matrix, values_array, ia, ja, x, e, res)

        write(*,*) '----------- Complex FEAST test END: -----------'
        write(*,*) " "

    end subroutine sparse_zfeast_test

    subroutine sparse_dfeast_test()
        
        implicit none
        double precision, allocatable :: test_matrix(:,:)
        CHARACTER*1 :: jobz, uplo            
        integer :: i, n, lda, lwork, info
        double precision, allocatable :: eigenvalues(:), rwork(:)
        double precision, allocatable ::  work(:,:)
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
        
        write(*,*) '----------- Double precision FEAST test START small matrix: -----------'
        write(*,*) " "

        n = 5
        allocate (test_matrix(n,n))
        test_matrix = 0.0d0
        test_matrix(1,1) = 1.0d0
        test_matrix(2,2) = 0.0d0
        test_matrix(3,3) = 2.0d0
        test_matrix(4,4) = 0.0d0
        test_matrix(5,5) = 3.0d0
        test_matrix(1,4) = 5.0d0
        test_matrix(2,3) = 6.0d0
        test_matrix(3,5) = 7.0d0
        
        !standard diagonalziation
        jobz='V' !N - only eigenvalues, V - eigenvalues + eigenvectors
        uplo='U' !upper hermitean matrix part passed to routine
        n=5
        lda=n
        lwork=3*n-1  
        allocate(work(lwork,lwork),rwork(3*n-2),eigenvalues(n))
        !call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
        call dsyev(jobz, uplo, n, test_matrix, lda, eigenvalues, work, lwork, info)
        if (info /= 0) then
            write(*,*) 'problem with dsyev, info=', info
        end if 
            
        do i=1,n
            write(*,*) i, eigenvalues(i) 
        enddo     
        deallocate(work,rwork,eigenvalues)

        write(*,*) 'dsyev diagonalization finished, FEAST start'

        !FEAST diagonalization
        ! we need to store diagonal zeros as well !!!
        allocate( values_array(8), ia(6), ja(8) )
        
        values_array(1) = 1.0d0
        values_array(2) = 5.0d0
        values_array(3) = 0.0d0
        values_array(4) = 6.0d0
        values_array(5) = 2.0d0
        values_array(6) = 7.0d0
        values_array(7) = 0.0d0
        values_array(8) = 3.0d0
        
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
        !write(*,*) 'Windows 11 new feature: feast might work only for Relase, not Debug!'
        write(*,*) 'Before dfeast_scsrev... '
        call dfeast_scsrev(uplo, n, values_array, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
        write(*,*) 'eps_out= ', epsout
        write(*,*) 'loop= ', loop
        write(*,*) 'dfeast_scsrev info=', info
        write(*,*) 'After dfeast_scsrev... '
            
        if (info /= 0) then
        write(*,*) 'problem with dfeast_scsrev, info=', info
        end if 
        
        write(*,*) 'no of eigenvalues found= ', m
        do i = 1, m
            write(*,*) i, e(i)
        end do
                    
        deallocate(test_matrix, values_array, ia, ja, x, e, res)

        write(*,*) '----------- Double precision FEAST test END -----------'
        write(*,*) " "

    end subroutine sparse_dfeast_test


end module tests_module