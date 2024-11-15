include "mkl_vsl.f90" !necessary for intel mkl random numbers

 
module omp_mkl_test
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

    end subroutine omp_mkl_small_test

end module omp_mkl_test

module spin_systems
    contains
    subroutine H_XXX_full()
        use omp_lib
        use mkl_vsl
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
end module spin_systems

program spin_code
    use spin_systems
    implicit none

    call H_XXX_full()

end program spin_code